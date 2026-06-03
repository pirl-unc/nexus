"""
The purpose of this python3 script is to calculate the Cancer Cell Fraction (CCF)
for somatic variants given an aligned tumor BAM file, a list of variants in Variant
Grammar format, a copy number segmentation TSV, and tumor purity.

CCF equation
------------
    CCF = VAF * [purity * CN_tumor + (1 - purity) * CN_normal] / (purity * multiplicity)

where:
    VAF          = alt / total                          (per-variant, see below)
    CN_tumor     = total copy number at the locus       (from the copy-number TSV)
    CN_normal    = germline copy number at the locus    (autosomes: --normal-ploidy,
                                                         default 2; sex chromosomes
                                                         depend on --sex, see below)
    purity       = tumor purity                         (user-supplied, in (0, 1])
    multiplicity = number of tumor chromosomes carrying
                   the variant                          (estimated per-variant; see below)

Multiplicity estimation (Tarabichi et al. 2021, Nature Methods, Box 1)
-------------------------------------------------------------------------------
Multiplicity is determined by the data once VAF, purity, and local copy number are known.
The paper's formula is:

    m_raw = (VAF / purity) * [purity * CN_tumor + (1 - purity) * CN_normal]
    m     = max(1, round(m_raw))
"""


import argparse
import math
import pandas as pd
import pysam
import sys
from pathlib import Path
from typing import Tuple


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate cancer cell fraction (CCF) for somatic DNA variants."
    )
    parser.add_argument(
        "--bam-file",
        dest="bam_file",
        type=Path,
        required=True,
        help="Tumor BAM file.",
    )
    parser.add_argument(
        "--variants-tsv-file",
        dest="variants_tsv_file",
        type=Path,
        required=True,
        help="Somatic DNA variants TSV file. Expected columns: variant_id, "
             "chromosome_1, position_1, strand_1, operation_1, chromosome_2, position_2, strand_2, operation_2, "
             "sequence, position_1_read_count_alternate_allele, position_2_read_count_alternate_allele."
    )
    parser.add_argument(
        "--copy-number-tsv-file",
        dest="copy_number_tsv_file",
        type=Path,
        required=True,
        help="Copy number TSV file. Expected columns: chromosome, start, end, "
             "copy_number, major_copy_number, minor_copy_number."
    )
    parser.add_argument(
        "--tumor-purity",
        dest="tumor_purity",
        type=float,
        required=True,
        help="Tumor purity in (0, 1]."
    )
    parser.add_argument(
        "--sex",
        dest="sex",
        type=str,
        required=True,
        choices=["male", "female"],
        help="Sex. Sets the germline (normal-cell) copy number on the sex chromosomes "
             "(chrX: 1 (male) / 2 (female); chrY: 1 (male) / 0 (female)). This information is used for CN_normal "
             "in the CCF and multiplicity formulas and the off-segment fallback genotype. Autosomes use --normal-ploidy."
    )
    parser.add_argument(
        "--normal-ploidy",
        dest="normal_ploidy",
        type=int,
        default=2,
        help="Normal (germline) copy number on the autosomes (default: 2). Sex chromosomes are set by --sex, not this value.",
    )
    parser.add_argument(
        "--min-base-quality",
        dest="min_base_quality",
        type=int,
        default=0,
        help="Minimum base quality for pileup depth at substitution loci (default: 0)."
    )
    parser.add_argument(
        "--min-mapping-quality",
        dest="min_mapping_quality",
        type=int,
        default=0,
        help="Minimum mapping quality for reads counted in depth (default: 0)."
    )
    parser.add_argument(
        "--anchor-bp",
        dest="anchor_bp",
        type=int,
        default=20,
        help="Flanking matched bases required on each side of a breakpoint for a read to count as "
             "reference-spanning at junction variants (default: 20)."
    )
    parser.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=4,
        help="Number of htslib threads for BGZF decompression of the phased BAM (default: 4)."
    )
    parser.add_argument(
        "--output-tsv-file",
        dest="output_tsv_file",
        type=Path,
        required=True,
        help="Output TSV file."
    )
    return parser.parse_args()


def _read_passes(read, min_mapping_quality: int) -> bool:
    return (
        not read.is_secondary
        and not read.is_supplementary
        and not read.is_duplicate
        and not read.is_unmapped
        and read.mapping_quality >= min_mapping_quality
    )


def pileup_depth(
        bam: pysam.AlignmentFile,
        chromosome: str,
        position: int,
        min_base_quality: int,
        min_mapping_quality: int
) -> int:
    """
    Fetch the number of reads covering 1-based `position` (reference + alternate), honoring quality filters.
    """
    try:
        coverage = bam.count_coverage(
            contig=chromosome,
            start=position - 1,
            stop=position,
            quality_threshold=min_base_quality,
            read_callback=lambda r: _read_passes(r, min_mapping_quality)
        )
    except (ValueError, KeyError):
        return 0
    # count_coverage returns four arrays (A, C, G, T) of length stop-start
    return int(sum(arr[0] for arr in coverage))


def reference_spanning_depth(
        bam: pysam.AlignmentFile,
        chromosome: str,
        position: int,
        anchor_bp: int,
        min_mapping_quality: int
) -> int:
    """
    Fetches the number of reads whose alignment continuously spans 1-based `position` with `anchor_bp` flanking match.

    A read is counted only if a single gapless aligned block covers [position - anchor, position + anchor].
    This excludes reads carrying a deletion as a CIGAR `D` (the block splits at the gap) and
    reads soft-clipped at the junction (the block ends at the breakpoint) -- i.e. it counts
    reference-supporting reads, not alt-supporting ones. Measured per breakpoint, so it works
    regardless of the event's size.
    """
    window_start = (position - 1) - anchor_bp        # 0-based inclusive
    window_end = (position - 1) + anchor_bp + 1      # 0-based exclusive
    count = 0
    try:
        for read in bam.fetch(chromosome, max(0, window_start), window_end):
            if not _read_passes(read, min_mapping_quality):
                continue
            for block_start, block_end in read.get_blocks(): # 0-based, end-exclusive
                if block_start <= window_start and block_end >= window_end:
                    count += 1
                    break
    except (ValueError, KeyError):
        return 0
    return count


def calculate_normal_copy_number(chromosome: str, sex: str, normal_ploidy: int) -> float:
    if sex == 'male':
        if chromosome in ['chrX', 'chrY', 'X', 'Y']:
            return 1.0
        return float(normal_ploidy)
    if sex == 'female':
        assert chromosome not in ['chrY', 'Y']
        return float(normal_ploidy)
    raise Exception("Unexpected to reach here.")


def fetch_local_copy_numbers(
        df_copy_numbers: pd.DataFrame,
        chromosome: str,
        position: int,
        sex: str,
        normal_ploidy: int
) -> Tuple[float, float, float]:
    df_matched = df_copy_numbers.loc[
        (df_copy_numbers['chromosome'] == chromosome) &
        (df_copy_numbers['start'] <= position) &
        (df_copy_numbers['end'] >= position),
        :
    ]
    if len(df_matched) == 0:
        normal_cn = calculate_normal_copy_number(
            chromosome=chromosome,
            sex=sex,
            normal_ploidy=normal_ploidy
        )
        return (normal_cn, math.ceil(normal_cn / 2.0), math.floor(normal_cn / 2.0))
    else:
        assert len(df_matched) == 1
        return (float(df_matched['copy_number'].values[0]),
                float(df_matched['major_copy_number'].values[0]),
                float(df_matched['minor_copy_number'].values[0]))


def calculate_vaf(
        df_variants: pd.DataFrame,
        bam_file: Path,
        min_base_quality: int,
        min_mapping_quality: int,
        anchor_bp: int,
        num_threads: int
) -> pd.DataFrame:
    vaf_values = []
    num_total_reads_values = []
    num_ref_reads_values = []
    with pysam.AlignmentFile(str(bam_file), "rb", threads=num_threads) as bam:
        for _, row in df_variants.iterrows():
            variant_type = row['variant_type']
            chromosome_1 = row['chromosome_1']
            position_1 = row['position_1']
            chromosome_2 = row['chromosome_2']
            position_2 = row['position_2']
            num_reads_1 = row['position_1_read_count_alternate_allele']
            num_reads_2 = row['position_2_read_count_alternate_allele']

            assert num_reads_1 > 0 and num_reads_2 > 0

            if num_reads_1 >= num_reads_2:
                num_alt_reads = num_reads_1
            else:
                num_alt_reads = num_reads_2

            if variant_type in ['SNV', 'MNV', 'INS']:
                d1 = pileup_depth(
                    bam=bam,
                    chromosome=chromosome_1,
                    position=position_1,
                    min_base_quality=min_base_quality,
                    min_mapping_quality=min_mapping_quality
                )
                d2 = pileup_depth(
                    bam=bam,
                    chromosome=chromosome_2,
                    position=position_2,
                    min_base_quality=min_base_quality,
                    min_mapping_quality=min_mapping_quality
                )
                num_total_reads = max((d1 + d2) / 2.0, num_alt_reads)
                num_reference_reads = num_total_reads - num_alt_reads
            else:
                d1 = reference_spanning_depth(
                    bam=bam,
                    chromosome=chromosome_1,
                    position=position_1,
                    anchor_bp=anchor_bp,
                    min_mapping_quality=min_mapping_quality
                )
                d2 = reference_spanning_depth(
                    bam=bam,
                    chromosome=chromosome_2,
                    position=position_2,
                    anchor_bp=anchor_bp,
                    min_mapping_quality=min_mapping_quality
                )
                num_reference_reads = (d1 + d2) / 2.0
                num_total_reads = num_alt_reads + num_reference_reads

            vaf = num_alt_reads / num_total_reads if num_total_reads > 0 else float("nan")

            vaf_values.append(vaf)
            num_total_reads_values.append(num_total_reads)
            num_ref_reads_values.append(num_reference_reads)

    df_variants['vaf'] = vaf_values
    df_variants['num_total_reads'] = num_total_reads_values
    df_variants['num_ref_reads'] = num_ref_reads_values

    return df_variants


def standardize_variant(
        chromosome_1: str,
        position_1: int,
        operation_1: str,
        chromosome_2: str,
        position_2: int,
        operation_2: str
) -> Tuple[str, int, str, str, int, str]:
    if chromosome_1 != chromosome_2:
        return chromosome_1, position_1, operation_1, chromosome_2, position_2, operation_2
    else:
        if position_1 < position_2:
            return chromosome_1, position_1, operation_1, chromosome_2, position_2, operation_2
        else:
            return chromosome_2, position_2, operation_2, chromosome_1, position_1, operation_1


def load_dna_variants(tsv_file: Path) -> pd.DataFrame:
    df_variants = pd.read_csv(tsv_file, sep="\t")
    required_cols = [
        "variant_id",
        "chromosome_1",
        "position_1",
        "strand_1",
        "operation_1",
        "chromosome_2",
        "position_2",
        "strand_2",
        "operation_2",
        "sequence",
        "position_1_read_count_alternate_allele",
        "position_2_read_count_alternate_allele"
    ]

    missing = [c for c in required_cols if c not in df_variants.columns]
    if missing:
        raise Exception(f"Missing required columns in --variants-tsv-file: {missing}")

    # Infer the variant types
    variant_types = []
    for _, row in df_variants.iterrows():
        chromosome_1 = str(row["chromosome_1"])
        position_1 = int(row["position_1"])
        operation_1 = str(row["operation_1"])
        chromosome_2 = str(row["chromosome_2"])
        position_2 = int(row["position_2"])
        operation_2 = str(row["operation_2"])

        # Standardize the variant
        (chromosome_1,
         position_1,
         operation_1,
         chromosome_2,
         position_2,
         operation_2) = standardize_variant(
            chromosome_1=chromosome_1,
            position_1=position_1,
            operation_1=operation_1,
            chromosome_2=chromosome_2,
            position_2=position_2,
            operation_2=operation_2
        )

        sequence = str(row["sequence"])

        # Determine the variant type
        variant_type = ''
        if chromosome_1 != chromosome_2:
            variant_type = 'TRA'
        else:
            if operation_1 == 'D' and operation_2 == 'U':
                if len(sequence) == 1 and abs(position_2 - position_1) == 2:
                    variant_type = 'SNV'
                if len(sequence) >= 2 and len(sequence) == abs(position_2 - position_1) - 1:
                    variant_type = 'MNV'
                if len(sequence) > 0  and abs(position_2 - position_1) == 1:
                    variant_type = 'INS'
                if len(sequence) == 0 and abs(position_2 - position_1) >= 2:
                    variant_type = 'DEL'
            if operation_1 == 'D' and operation_2 == 'D':
                variant_type = 'INV' # head-to-head
            if operation_1 == 'U' and operation_2 == 'U':
                variant_type = 'INV' # tail-to-tail
            if operation_1 == 'U' and operation_2 == 'D':
                variant_type = 'DUP'

        if variant_type == '':
            raise ValueError(
                f"Could not infer variant type for variant_id={row['variant_id']} "
                f"({chromosome_1}:{position_1}:{operation_1} / {chromosome_2}:{position_2}:{operation_2}, "
                f"seq_len={len(sequence)})"
            )

        variant_types.append(variant_type)

    df_variants["variant_type"] = variant_types

    return df_variants


def load_copy_numbers(tsv_file: Path) -> pd.DataFrame:
    df_copy_number = pd.read_csv(tsv_file, sep="\t")
    for c in ("chromosome", "start", "end", "copy_number", "major_copy_number", "minor_copy_number"):
        if c not in df_copy_number.columns:
            raise Exception(f"Missing required column in --copy-number-tsv-file: {c}")
    return df_copy_number


def run():
    args = parse_args()

    # Step 1. Check inputs
    if not (0 < args.tumor_purity <= 1):
        sys.exit(f"--tumor-purity must be in (0, 1]; got {args.tumor_purity}")

    # Step 2. Load somatic DNA variants data
    df_variants = load_dna_variants(tsv_file=args.variants_tsv_file)

    # Step 3. Load copy number data
    df_copy_number = load_copy_numbers(tsv_file=args.copy_number_tsv_file)

    # Step 4. Calculate variant allele fraction (VAF) for each variant
    df_variants = calculate_vaf(
        df_variants=df_variants,
        bam_file=args.bam_file,
        min_base_quality=args.min_base_quality,
        min_mapping_quality=args.min_mapping_quality,
        anchor_bp=args.anchor_bp,
        num_threads=args.num_threads
    )

    # Step 5. Calculate cancer cell fraction (CCF) for each variant
    tumor_purity = args.tumor_purity
    ccf_values = []
    mutation_multiplicity_values = []
    tumor_cn_values = []
    major_cn_values = []
    minor_cn_values = []
    normal_cn_values = []
    for _, row in df_variants.iterrows():
        chromosome_1 = str(row["chromosome_1"])
        position_1 = int(row["position_1"])
        chromosome_2 = str(row["chromosome_2"])
        position_2 = int(row["position_2"])
        vaf = float(row['vaf'])

        # Get the local copy numbers
        cn_1, major_cn_1, minor_cn_1 = fetch_local_copy_numbers(
            df_copy_numbers=df_copy_number,
            chromosome=chromosome_1,
            position=position_1,
            sex=args.sex,
            normal_ploidy=args.normal_ploidy
        )
        cn_2, major_cn_2, minor_cn_2 = fetch_local_copy_numbers(
            df_copy_numbers=df_copy_number,
            chromosome=chromosome_2,
            position=position_2,
            sex=args.sex,
            normal_ploidy=args.normal_ploidy
        )
        tumor_cn = (float(cn_1) + float(cn_2)) / 2.0
        major_cn = (float(major_cn_1) + float(major_cn_2)) / 2.0
        minor_cn = (float(minor_cn_1) + float(minor_cn_2)) / 2.0

        # Calculate the normal ploidy
        normal_cn_1 = calculate_normal_copy_number(
            chromosome=chromosome_1,
            sex=args.sex,
            normal_ploidy=args.normal_ploidy
        )
        normal_cn_2 = calculate_normal_copy_number(
            chromosome=chromosome_2,
            sex=args.sex,
            normal_ploidy=args.normal_ploidy
        )
        if normal_cn_1 != normal_cn_2:
            normal_cn = (normal_cn_1 + normal_cn_2) / 2.0
        else:
            normal_cn = normal_cn_1

        # Calculate m multiplicity
        m_raw = (vaf / tumor_purity) * ((tumor_purity * tumor_cn) + (normal_cn * (1 - tumor_purity)))
        m = max(1, math.floor(m_raw + 0.5))

        # Calculate CCF
        ccf = (vaf / (m * tumor_purity)) * ((tumor_purity * tumor_cn) + (normal_cn * (1 - tumor_purity)))

        ccf_values.append(ccf)
        mutation_multiplicity_values.append(m)
        tumor_cn_values.append(tumor_cn)
        major_cn_values.append(major_cn)
        minor_cn_values.append(minor_cn)
        normal_cn_values.append(normal_cn)

    df_variants['ccf'] = ccf_values
    df_variants["mutation_multiplicity"] = mutation_multiplicity_values
    df_variants['tumor_copy_number'] = tumor_cn_values
    df_variants['tumor_major_copy_number'] = major_cn_values
    df_variants['tumor_minor_copy_number'] = minor_cn_values
    df_variants['normal_copy_number'] = normal_cn_values

    df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)

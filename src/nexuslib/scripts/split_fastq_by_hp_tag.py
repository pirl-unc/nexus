"""
Split a fastq.gz file into three fastq.gz files (haplotype 1, haplotype 2, and unknown)
based on the HP tag in a phased BAM file.

Rules:
  - HP=1 only across all alignments of a read -> hap1 output
  - HP=2 only across all alignments of a read -> hap2 output
  - No HP tag on any alignment of a read      -> unknown output
  - Conflicting HP tags across primary +
    supplementary/secondary alignments         -> unknown output
  - Read present in fastq but absent from BAM  -> controlled by
    --missing-from-bam (default: exclude)

Usage:
    python split_fastq_by_hp_tag.py \
        --bam-file phased.bam \
        --fastq-file reads.fastq.gz \
        --out-hap1-fastq-file reads.hap1.fastq.gz \
        --out-hap2-fastq-file reads.hap2.fastq.gz \
        --out-hap-unknown-fastq-file reads.unknown.fastq.gz \
        [--missing-from-bam {exclude,unknown}]
"""


import argparse
import pysam
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set, Tuple


def parse_args():
    parser = argparse.ArgumentParser(
        description="Split a fastq.gz file by HP tag from a phased BAM file."
    )
    parser.add_argument(
        "--bam-file",
        required=True,
        type=Path,
        dest="bam_file",
        help="Phased BAM file with HP tags."
    )
    parser.add_argument(
        "--fastq-file",
        required=True,
        type=Path,
        dest="fastq_file",
        help="Original fastq.gz file."
    )
    parser.add_argument(
        "--out-hap1-fastq-file",
        required=True,
        type=Path,
        dest="out_hap1_fastq_file",
        help="Output fastq.gz for haplotype 1."
    )
    parser.add_argument(
        "--out-hap2-fastq-file",
        required=True,
        type=Path,
        dest="out_hap2_fastq_file",
        help="Output fastq.gz for haplotype 2."
    )
    parser.add_argument(
        "--out-hap-unknown-fastq-file",
        required=True,
        type=Path,
        dest="out_hapunknown_fastq_file",
        help=(
            "Output fastq.gz for reads with no HP tag on any alignment, or with "
            "conflicting HP tags across alignments."
        )
    )
    parser.add_argument(
        "--missing-from-bam",
        dest="missing_from_bam",
        choices=["exclude", "unknown"],
        default="exclude",
        help=(
            "How to handle reads in the fastq that are absent from the BAM. "
            "'exclude' (default) drops them. 'unknown' routes them to --out-hap-unknown."
        )
    )
    parser.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=4,
        help="Number of htslib threads for BGZF decompression of the phased BAM (default: 4)."
    )
    return parser.parse_args()


def _fastq_record_bytes(entry):
    """
    Serialize a pysam.FastxFile entry back to a 4-line FASTQ record (as bytes,
    ready to write to a pysam.BGZFile). The original header is preserved: the
    bare read id (``entry.name``) plus its comment, if any.
    """
    header = entry.name if entry.comment is None else "%s %s" % (entry.name, entry.comment)
    return ("@%s\n%s\n+\n%s\n" % (header, entry.sequence, entry.quality)).encode()


def collect_hp_tags(
        bam_file: str,
        num_threads: int = 1
) -> Tuple[Set[str], Dict[str, Set[int]]]:
    """
    Walk every alignment in the BAM (primary, supplementary, secondary, unmapped)
    and aggregate the set of HP tag values seen per read name.

    Parameters:
        bam_file    :   BAM file.
        num_threads :   Number of threads.

    Returns:
        seen_reads  :   Set[str] of all read names observed in the BAM file.
        hp_tags     :   Dict[str, Set[int]] where the key is read name and the value is Set of HP tag values.
                        If a key is missing or if the set is empty => no HP tag on any alignment.
    """
    seen_reads = set()
    hp_tags = defaultdict(set)
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False, threads=num_threads) as bam:
        for aln in bam:
            seen_reads.add(aln.query_name)
            try:
                hp_tags[aln.query_name].add(int(aln.get_tag("HP")))
            except KeyError:
                # No HP tag on this alignment
                pass
    return seen_reads, hp_tags


def bucket_for(
        read_name: str,
        seen_reads: Set[str],
        hp_tags: Dict[str, Set[int]],
        missing_from_bam: str
) -> str:
    """Return one of: 'hap1', 'hap2', 'unknown', 'exclude'."""
    if read_name not in seen_reads:
        # Read not in BAM at all -> behavior controlled by --missing-from-bam
        return missing_from_bam  # 'exclude' or 'unknown'
    tags = hp_tags.get(read_name, set())
    if tags == {1}:
        return "hap1"
    if tags == {2}:
        return "hap2"
    # Empty (no HP) OR conflicting (e.g., {1, 2}) -> unknown
    return "unknown"


def run():
    args = parse_args()

    # Step 1. Collect HP tags per read name from BAM
    print("Reading HP tags from BAM: %s (threads=%i)" % (args.bam_file, args.num_threads))
    seen_reads, hp_tags = collect_hp_tags(
        bam_file=args.bam_file,
        num_threads=args.num_threads
    )
    print("\ti unique read names in BAM file." % len(seen_reads))
    print("\t%i read names with at least one HP tag." % len(hp_tags))

    # Step 2. Make output directories
    args.out_hap1_fastq_file.parent.mkdir(parents=True, exist_ok=True)
    args.out_hap2_fastq_file.parent.mkdir(parents=True, exist_ok=True)
    args.out_hapunknown_fastq_file.parent.mkdir(parents=True, exist_ok=True)

    # Step 3. Read the input fastq with pysam and write to the three outputs.
    print("Splitting fastq: %s" % args.fastq_file)
    print("\t--missing-from-bam = %s" % args.missing_from_bam)
    n_hap1 = n_hap2 = 0
    n_unknown_conflict_or_missing_hp = 0
    n_unknown_not_in_bam = n_excluded_not_in_bam = 0
    with pysam.FastxFile(args.fastq_file) as fin, \
         pysam.BGZFile(args.out_hap1_fastq_file, "wb") as out1, \
         pysam.BGZFile(args.out_hap2_fastq_file, "wb") as out2, \
         pysam.BGZFile(args.out_hapunknown_fastq_file, "wb") as out_unknown:
        for entry in fin:
            read_name = str(entry.name)
            b = bucket_for(
                read_name=read_name,
                seen_reads=seen_reads,
                hp_tags=hp_tags,
                missing_from_bam=args.missing_from_bam
            )
            if b == "hap1":
                out1.write(_fastq_record_bytes(entry))
                n_hap1 += 1
            elif b == "hap2":
                out2.write(_fastq_record_bytes(entry))
                n_hap2 += 1
            elif b == "unknown":
                out_unknown.write(_fastq_record_bytes(entry))
                if read_name not in seen_reads:
                    n_unknown_not_in_bam += 1
                else:
                    n_unknown_conflict_or_missing_hp += 1
            else:
                # 'exclude'
                n_excluded_not_in_bam += 1

    # Step 4. Report
    total = (
        n_hap1
        + n_hap2
        + n_unknown_conflict_or_missing_hp
        + n_unknown_not_in_bam
        + n_excluded_not_in_bam
    )
    print("Done.")
    print("\thap1-only reads                    = %i" % n_hap1)
    print("\thap2-only reads                    = %i" % n_hap2)
    print("\tunknown (no HP / conflicting HP)   = %i" % n_unknown_conflict_or_missing_hp)
    print("\tunknown (read absent from BAM)     = %i" % n_unknown_not_in_bam)
    print("\texcluded (read absent from BAM)    = %i" % n_excluded_not_in_bam)
    print("\ttotal reads processed              = %i" % total)

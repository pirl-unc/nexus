"""
Merge per-haplotype RNA assembly outputs (H1, H2, Hunknown) into:

  1. A single FASTA whose sequence IDs are prefixed with the haplotype tag
     (`h1_`, `h2_`, `hunknown_`) so downstream tools can distinguish the
     three contributing assemblies without losing the original IDs.
  2. A TSV indexed by the prefixed ID, carrying the sequence and the
     comma-joined list of supporting read names (joined from the per-
     haplotype reads.tsv on the original `transcript_id`).
  3. A gzipped FASTQ mirroring (1) — same prefixed IDs, same sequences,
     with a fixed per-base quality (default Phred 60) so downstream tools
     that want FASTQ input (e.g. aligners) can consume the merged assembly
     directly. Assembled contigs do not carry meaningful base qualities,
     so a constant value is the honest choice; --base-quality lets you
     tune it for tools that filter on Q.

The input reads.tsv schema (one row per read-to-transcript alignment):
  read_name,
  read_length, read_start, read_end, strand,
  transcript_id, transcript_length, transcript_start, transcript_end,
  num_residue_matches, frac_residue_matches, alignment_block_length,
  mapping_quality

Usage (typically invoked from a sibling bash wrapper per sample):
  python merge_scrna_assembly_files.py \\
      --h1-fasta        <h1.fa>       --h1-reads-tsv       <h1.reads.tsv>      \\
      --h2-fasta        <h2.fa>       --h2-reads-tsv       <h2.reads.tsv>      \\
      --hunknown-fasta  <hu.fa>       --hunknown-reads-tsv <hu.reads.tsv>      \\
      --output-fasta    <merged.fa>   --output-tsv         <merged.tsv>        \\
      --output-fastq-gz <merged.fq.gz> [--base-quality 60]
"""


import argparse
import gzip
import pandas as pd
import pysam
import sys
from pathlib import Path
from typing import Dict, List


HAPLOTYPE_PREFIX: Dict[str, str] = {
    "h1":       "h1_",
    "h2":       "h2_",
    "hunknown": "hunknown_",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Merge hap 1 / hap 2 / unknown hap RNA assembly FASTAs and their read support information "
            "into one prefixed FASTA file and one read support TSV file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        "--h1-fasta-file",
        required=True,
        type=Path,
        dest="h1_fasta",
        help="Haplotype 1 RNA assembly FASTA file"
    )
    p.add_argument(
        "--h2-fasta-file",
        required=True,
        type=Path,
        dest="h2_fasta",
        help="Haplotype 2 RNA assembly FASTA file"
    )
    p.add_argument(
        "--hunknown-fasta-file",
        required=True,
        type=Path,
        dest="hunknown_fasta",
        help="Unknown haplotype RNA assembly FASTA file"
    )
    p.add_argument(
        "--h1-reads-tsv-file",
        required=True,
        type=Path,
        dest="h1_reads_tsv",
        help="Haplotype 1 RNA assembly read support file. Expected columns: 'read_name', 'transcript_id'. "
             "A transcript ID with multiple read support should have multiple rows."
    )
    p.add_argument(
        "--h2-reads-tsv-file",
        required=True,
        type=Path,
        dest="h2_reads_tsv",
        help="Haplotype 2 RNA assembly read support file. Expected columns: 'read_name', 'transcript_id'. "
             "A transcript ID with multiple read support should have multiple rows."
    )
    p.add_argument(
        "--hunknown-reads-tsv-file",
        required=True,
        type=Path,
        dest="hunknown_reads_tsv",
        help="Unknown haplotype RNA assembly read support file. Expected columns: 'read_name', 'transcript_id'. "
             "A transcript ID with multiple read support should have multiple rows."
    )
    p.add_argument(
        "--output-fasta-file",
        required=True,
        type=Path,
        dest="output_fasta_file",
        help="Output FASTA file."
    )
    p.add_argument(
        "--output-tsv-file",
        required=True,
        type=Path,
        dest="output_tsv_file",
        help="Output TSV file."
    )
    p.add_argument(
        "--output-fastq-file",
        required=True,
        type=Path,
        dest="output_fastq_file",
        help="Output FASTQ.GZ that mirrors --output-fasta-file with a constant per-base quality (see --base-quality)."
    )
    p.add_argument(
        "--base-quality",
        type=int,
        default=60,
        dest="base_quality",
        help="Phred quality assigned to every base in the FASTQ output. Valid range: 0-93 (Phred+33 encoding)."
    )
    return p.parse_args()


def load_reads_index(reads_tsv_file: Path) -> Dict[str, List[str]]:
    """
    Load reads.

    Parameters:
        reads_tsv_file  :   Reads TSV file path.

    Returns:
        reads_dict      :   Dict[transcript ID, List[read IDs]]
    """
    df_reads = pd.read_csv(reads_tsv_file, sep="\t")
    read_dict: Dict[str, List[str]] = {}
    for tid, sub in df_reads.groupby("transcript_id", sort=False):
        read_dict[str(tid)] = sub["read_name"].astype(str).tolist()
    return read_dict


def run():
    args = parse_args()

    # Step 1. Check inputs
    # Validate base quality. Phred+33 is printable ASCII 33..126 (i.e. quality scores 0..93).
    # Out-of-range values would produce an invalid quality string.
    if not (0 <= args.base_quality <= 93):
        sys.exit(
            f"--base-quality must be in [0, 93] (Phred+33); "
            f"got {args.base_quality}"
        )
    quality_char = chr(args.base_quality + 33)

    # Step 2. Make output directories
    args.output_fasta_file.parent.mkdir(parents=True, exist_ok=True)
    args.output_tsv_file.parent.mkdir(parents=True, exist_ok=True)
    args.output_fastq_file.parent.mkdir(parents=True, exist_ok=True)

    # Step 3. Write the merged output FASTA and FASTQ files
    inputs = [
        ("h1", args.h1_fasta, args.h1_reads_tsv),
        ("h2", args.h2_fasta, args.h2_reads_tsv),
        ("hunknown", args.hunknown_fasta, args.hunknown_reads_tsv)
    ]
    data = {
        'assembled_transcript_name': [],
        'sequence': [],
        'read_names': [],
        'num_read_names': []
    }
    seen_read_ids = set()
    with open(args.output_fasta_file, "w") as fasta_out_fh, gzip.open(args.output_fastq_file, "wt") as fastq_out_fh:
        for haplotype, fasta_file, reads_tsv_file in inputs:
            print("Merge:")
            print("\tHaplotype: %s" % haplotype)
            print("\tFASTA file: %s" % fasta_file)
            print("\tReads TSV file: %s" % reads_tsv_file)
            prefix = HAPLOTYPE_PREFIX[haplotype]
            reads_dict = load_reads_index(reads_tsv_file=reads_tsv_file)
            with pysam.FastxFile(str(fasta_file)) as fh:
                for record in fh:
                    transcript_id = str(record.name)
                    sequence = str(record.sequence)

                    new_transcript_id = "%s%s" % (prefix, transcript_id)

                    if new_transcript_id in seen_read_ids:
                        raise Exception("We have previously seen transcript ID: %s" % transcript_id)
                    seen_read_ids.add(new_transcript_id)

                    fasta_out_fh.write(">%s\n%s\n" % (new_transcript_id, sequence))
                    fastq_out_fh.write("@%s\n%s\n+\n%s\n" % (new_transcript_id, sequence, quality_char * len(sequence)))

                    if transcript_id not in reads_dict:
                        raise Exception("%s is expected to be in %s" % (transcript_id, reads_tsv_file))

                    read_names = reads_dict[transcript_id]

                    data['assembled_transcript_name'].append(new_transcript_id)
                    data['sequence'].append(sequence)
                    data['read_names'].append(";".join(read_names))
                    data['num_read_names'].append(len(read_names))

    df_reads = pd.DataFrame(data)
    df_reads.to_csv(args.output_tsv_file, sep="\t", index=False)

    print("Done.")
    print("\tTotal number of transcripts: %i" % len(df_reads))

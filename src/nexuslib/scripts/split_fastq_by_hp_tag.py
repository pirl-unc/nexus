"""
Split a fastq.gz file into two fastq.gz files (haplotype 1 and haplotype 2)
based on the HP tag in a phased BAM file.

Rules:
  - HP=1 only across all alignments of a read -> hap1 output
  - HP=2 only across all alignments of a read -> hap2 output
  - No HP tag on any alignment of a read     -> both outputs
  - Conflicting HP tags across primary +
    supplementary/secondary alignments        -> both outputs
  - Read present in fastq but absent from BAM -> controlled by
    --missing-from-bam (default: exclude)

Usage:
    python split_fastq_by_hp_tag.py \
        --bam phased.bam \
        --fastq reads.fastq.gz \
        --out-hap1 reads.hap1.fastq.gz \
        --out-hap2 reads.hap2.fastq.gz \
        [--missing-from-bam {exclude,both}]
"""


import argparse
import dnaio
import pysam
import subprocess
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Split a fastq.gz file by HP tag from a phased BAM file."
    )
    parser.add_argument("--bam-file", required=True, dest="bam", help="Phased BAM file with HP tags.")
    parser.add_argument("--fastq-file", required=True, dest="fastq_file", help="Original fastq.gz file.")
    parser.add_argument("--out-hap1", required=True, dest="out_hap_1", help="Output fastq.gz for haplotype 1.")
    parser.add_argument("--out-hap2", required=True, dest="out_hap_2", help="Output fastq.gz for haplotype 2.")
    parser.add_argument(
        "--missing-from-bam",
        choices=["exclude", "both"],
        default="exclude",
        help=(
            "How to handle reads in the fastq that are absent from the BAM. "
            "'exclude' (default) drops them. 'both' writes them to both outputs."
        ),
    )
    return parser.parse_args()


def collect_hp_tags(bam_file):
    """
    Walk every alignment in the BAM (primary, supplementary, secondary, unmapped)
    and aggregate the set of HP tag values seen per read name.

    Returns:
        seen_reads: set[str] of all read names observed in the BAM
        hp_tags:    dict[str, set[int]] of HP tag values per read name
                    (missing key OR empty set => no HP tag on any alignment)
    """
    seen_reads = set()
    hp_tags = defaultdict(set)
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
        for aln in bam:
            seen_reads.add(aln.query_name)
            try:
                hp_tags[aln.query_name].add(int(aln.get_tag("HP")))
            except KeyError:
                pass  # no HP tag on this alignment
    return seen_reads, hp_tags


def bucket_for(read_name, seen_reads, hp_tags, missing_from_bam):
    """Return one of: 'hap1', 'hap2', 'both', 'exclude'."""
    if read_name not in seen_reads:
        # Read not in BAM at all -> behavior controlled by --missing-from-bam
        return missing_from_bam  # 'exclude' or 'both'
    tags = hp_tags.get(read_name, set())
    if tags == {1}:
        return "hap1"
    if tags == {2}:
        return "hap2"
    # Empty (no HP) OR conflicting (e.g., {1, 2}) -> both
    return "both"


def run():
    args = parse_args()

    # Step 1. Collect HP tags per read name from BAM.
    print("Reading HP tags from BAM: %s" % args.bam)
    seen_reads, hp_tags = collect_hp_tags(args.bam)
    print("  %d unique read names in BAM" % len(seen_reads))
    print("  %d read names with at least one HP tag" % len(hp_tags))

    # Step 2. Stream the input fastq via pigz and write to the two outputs.
    print("Splitting fastq: %s" % args.fastq_file)
    print("  --missing-from-bam = %s" % args.missing_from_bam)
    n_hap1 = n_hap2 = 0
    n_both_conflict_or_missing_hp = 0
    n_both_not_in_bam = n_excluded_not_in_bam = 0
    proc = subprocess.Popen(["pigz", "-dc", args.fastq_file], stdout=subprocess.PIPE)
    with dnaio.open(proc.stdout, fileformat="fastq") as fin, \
         dnaio.open(args.out_hap_1, mode="w") as out1, \
         dnaio.open(args.out_hap_2, mode="w") as out2:
        for r in fin:
            read_name = r.name.split()[0]
            b = bucket_for(read_name, seen_reads, hp_tags, args.missing_from_bam)
            if b == "hap1":
                out1.write(r)
                n_hap1 += 1
            elif b == "hap2":
                out2.write(r)
                n_hap2 += 1
            elif b == "both":
                out1.write(r)
                out2.write(r)
                if read_name not in seen_reads:
                    n_both_not_in_bam += 1
                else:
                    n_both_conflict_or_missing_hp += 1
            else:  # 'exclude'
                n_excluded_not_in_bam += 1
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError("pigz exited with code %d" % proc.returncode)

    # Step 3. Report.
    total = (
        n_hap1
        + n_hap2
        + n_both_conflict_or_missing_hp
        + n_both_not_in_bam
        + n_excluded_not_in_bam
    )
    print("Done.")
    print("\thap1-only reads                    = %d" % n_hap1)
    print("\thap2-only reads                    = %d" % n_hap2)
    print("\tboth (no HP / conflicting HP)      = %d" % n_both_conflict_or_missing_hp)
    print("\tboth (read absent from BAM)        = %d" % n_both_not_in_bam)
    print("\texcluded (read absent from BAM)    = %d" % n_excluded_not_in_bam)
    print("\ttotal reads processed              = %d" % total)

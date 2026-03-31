"""
The purpose of this python3 file is to convert an unaligned BAM file to a FASTQ.GZ file while filtering out
cells that do not meet a minimum molecule count.
"""


import argparse
import gzip
import pandas as pd
import pysam
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(description="Convert an unaligned BAM file to a FASTQ.GZ file while filtering for barcodes with sufficient reads.")
    parser.add_argument("--bam-file", dest='bam_file', type=str, help="BAM file. Make sure this is an unaligned BAM file.")
    parser.add_argument("--barcode-tag", dest='barcode_tag', type=str, help="Barcode tag (e.g. CB for PacBio single-cell RNA-seq BAM file).")
    parser.add_argument("--output-fastq-file", dest='output_fastq_file', type=str, help="Output FASTQ.GZ file.")
    parser.add_argument("--output-tsv-file", dest='output_tsv_file', type=str, help="Output TSV file.")
    parser.add_argument("--min-reads-per-barcode", dest='min_reads_per_barcode', type=int, help="Minimum number of reads per barcode. Determine this value by running nexus_find_scrna_barcode_knee.")
    parser.add_argument("--base-quality", dest='base_quality', type=int, default=60, help="Default base quality score to use when reads lack quality scores (default: 60).")
    parser.add_argument("--num-threads", dest='num_threads', type=int, default=4, help="Number of threads (default: 4).")
    arguments = parser.parse_args()
    return arguments


def run():
    args = parse_args()
    default_qual_char = chr(args.base_quality + 33)

    # Step 1. Count the number of reads for each barcode
    counter = Counter()
    with pysam.AlignmentFile(args.bam_file, "rb", check_sq=False, threads=args.num_threads) as bam:
        for read in bam:
            if read.has_tag(args.barcode_tag):
                cb = read.get_tag(args.barcode_tag)
                counter[cb] += 1

    # Step 2. Get barcodes with the minimum number of reads
    barcodes = set(b for b, c in counter.items() if c >= args.min_reads_per_barcode)

    # Step 3. Iterate through the eligible barcodes and output the reads to a single FASTQ.GZ file
    data = {
        'cell_barcode': [],
        'read_name': []
    }
    with pysam.AlignmentFile(args.bam_file, "rb", check_sq=False, threads=args.num_threads) as bam, gzip.open(args.output_fastq_file, "wt") as out:
        for read in bam:
            if read.has_tag(args.barcode_tag):
                cb = read.get_tag(args.barcode_tag)
                if cb in barcodes:
                    if read.query_qualities is not None:
                        qual_str = pysam.qualities_to_qualitystring(read.query_qualities)
                    else:
                        qual_str = default_qual_char * read.query_length
                    out.write(f"@{read.query_name}\n{read.query_sequence}\n+\n{qual_str}\n")
                    data['cell_barcode'].append(cb)
                    data['read_name'].append(read.query_name)
    df = pd.DataFrame(data)
    df.to_csv(args.output_tsv_file, sep='\t', index=False)

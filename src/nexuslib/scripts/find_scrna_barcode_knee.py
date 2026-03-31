"""
The purpose of this python3 script is to identify the knee point of an scRNA-seq barcode rank plot from an
unaligned BAM file. The resulting threshold can be used to set a minimum molecule count per cell, filtering
out low-count barcodes in downstream analysis.
"""


import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(description="Find the knee of the barcode rank plot from an unaligned scRNA-seq BAM file.")
    parser.add_argument("--bam-file", dest='bam_file', type=str, help="Input BAM file (unaligned).")
    parser.add_argument("--tag", dest='tag', type=str, default='CB', help="Cell barcode tag (default: 'CB').")
    parser.add_argument("--percentile", dest='percentile', type=float, default=99, help="Percentile method's percentile to determine the knee (default: 99').")
    parser.add_argument("--multiplier", dest='multiplier', type=float, default=0.1, help="Percentile method's multiplier to determine the knee (default: 0.1').")
    parser.add_argument("--output-tsv-file", dest='output_tsv_file', type=str, required=False, help="Output TSV file.")
    parser.add_argument("--output-png-file", dest='output_png_file', type=str, required=False, help="Output PNG file. If you specify this, the plot will not be shown via the terminal.")
    parser.add_argument("--num-threads", dest='num_threads', type=int, default=4, help="Number of threads (default: 4).")
    arguments = parser.parse_args()
    return arguments


def find_knee_percentile(barcodes, percentile=99, multiplier=0.1):
    counts = np.array(list(barcodes.values()))
    cutoff = np.percentile(counts, percentile) * multiplier
    passing = {bc: c for bc, c in barcodes.items() if c >= cutoff}
    return cutoff, passing


def run():
    # Step 1. Parse arguments
    args = parse_args()

    # Step 2. Count the cell barcodes
    barcodes = Counter() # key = cell barcode, value = molecule count
    with pysam.AlignmentFile(args.bam_file, "rb", check_sq=False, threads=args.num_threads) as bam:
        for read in bam:
            if read.has_tag(args.tag):
                cb = read.get_tag(args.tag)
                barcodes[cb] += 1
    print("%i unique barcodes in total." % len(barcodes))

    # Step 3. Save to TSV file (optional)
    if args.output_tsv_file is not None:
        print("Saving to TSV file.")
        df = pd.DataFrame.from_dict(barcodes, orient="index", columns=["count"]).reset_index()
        df.columns = ["barcode", "count"]
        df.to_csv(args.output_tsv_file, sep="\t", index=False)

    # Step 4. Identify the knee
    cutoff, passing = find_knee_percentile(
        barcodes=barcodes,
        percentile=args.percentile,
        multiplier=args.multiplier
    )
    print(f"Threshold: {cutoff:.0f} reads")
    print(f"{len(passing):,} barcodes passing")

    # Step 5. Plot knee plot
    vals = sorted(barcodes.values(), reverse=True)
    plt.plot(vals)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Barcode rank")
    plt.ylabel("Reads per barcode")
    plt.title(f"Barcode rank plot ({len(barcodes):,} total cells, {sum(barcodes.values()):,} total molecules)")

    if args.output_png_file is not None:
        plt.savefig(args.output_png_file, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()

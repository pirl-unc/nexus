"""
The purpose of this python3 file is to count cell and molecule support for each transcript from assembled from scRNA-seq.
"""


import argparse
import pandas as pd
import pysam
import re


def parse_args():
    parser = argparse.ArgumentParser(description="Count cell and molecule support for each transcript from assembled from scRNA-seq.")
    parser.add_argument("--fasta-file", dest='fasta_file', type=str, help="Transcriptome assembly FASTA file.")
    parser.add_argument("--bam-file", dest='bam_file', type=str, help="Unaligned BAM file.")
    parser.add_argument("--reads-tsv-file", dest='reads_tsv_file', type=str, help="Reads TSV file (expected columns: 'transcript_id', 'read_name').")
    parser.add_argument("--output-tsv-file", dest='output_tsv_file', type=str, help="Output TSV file.")
    parser.add_argument("--tag", dest='tag', type=str, default='CB', help="Cell barcode tag in the BAM file (default: 'CB').")
    arguments = parser.parse_args()
    return arguments


def process_transcript(transcript_model_id, molecule_ids, molecule_ids_to_cb):
    cell_barcodes = set()
    for mid in molecule_ids:
        mid_ = re.sub(r'_s\d+$', '', mid)
        if mid_ in molecule_ids_to_cb:
            cell_barcodes.add(molecule_ids_to_cb[mid_])
    return transcript_model_id, molecule_ids, cell_barcodes


def run():
    args = parse_args()

    # Step 1. Load the transcript model IDs
    with pysam.FastaFile(args.fasta_file) as fa:
        transcript_model_ids = set(fa.references)

    # Step 2. Load the molecule name / cell barcode pair information
    molecule_id_to_cb = {} # key = molecule ID, value = cell barcode
    cell_barcode_ids = set()
    with pysam.AlignmentFile(args.bam_file, "rb", check_sq=False) as bam:
        for read in bam:
            assert read.query_name not in molecule_id_to_cb, f"Duplicate read name: {read.query_name}"
            cell_barcode = read.get_tag(args.tag)
            molecule_id_to_cb[read.query_name] = cell_barcode
            cell_barcode_ids.add(cell_barcode)
    print(f"{len(molecule_id_to_cb)} molecule IDs identified in the raw BAM file.", flush=True)
    print(f"{len(cell_barcode_ids)} cell barcodes identified in the raw BAM file.", flush=True)

    # Step 3. Load the read support TSV and build the lookup
    df_reads = pd.read_csv(args.reads_tsv_file, sep="\t", dtype={"transcript_id": str})
    transcript_to_reads = df_reads.groupby("transcript_id")["read_name"].apply(set).to_dict() # read name is molecule ID

    # Step 4. Process transcripts
    print(f"Processing {len(transcript_model_ids)} transcripts...", flush=True)
    results = [
        process_transcript(tid, transcript_to_reads.get(tid, set()), molecule_id_to_cb)
        for tid in transcript_model_ids
    ]

    # Step 6. Assemble DataFrame, sort by num_cell_barcodes descending
    df_results = pd.DataFrame(
        [
            {
                "transcript_id": tid,
                "molecule_ids": ",".join(sorted(read_names)),
                "cell_barcodes": ",".join(sorted(cell_barcodes)),
                "num_molecules": len(read_names),
                "num_cells": len(cell_barcodes),
            }
            for tid, read_names, cell_barcodes in results
        ]
    )
    df_results = df_results.sort_values("num_cells", ascending=False).reset_index(drop=True)

    df_results.to_csv(args.output_tsv_file, sep="\t", index=False)

    print(f"Output written to {args.output_tsv_file}", flush=True)

"""
The purpose of this python3 file is to filter scRNA-seq assembly based on minimum cell and molecule support.
"""


import argparse
import pandas as pd
import pysam


def parse_args():
    parser = argparse.ArgumentParser(description="Filter scRNA-seq assembly based on minimum cell and molecule support.")
    parser.add_argument("--fasta-file", dest='fasta_file', type=str, help="Transcriptome assembly FASTA file.")
    parser.add_argument("--support-tsv-file", dest='support_tsv_file', type=str, help="Support TSV file (expected columns: 'transcript_id', 'num_molecules', 'num_cells').")
    parser.add_argument("--min-cells", dest='min_cells', type=int, help="Minimum cell count.")
    parser.add_argument("--min-molecules", dest='min_molecules', type=int, help="Minimum molecule count.")
    parser.add_argument("--output-fasta-file", dest='output_fasta_file', type=str, help="Output FASTA file.")
    parser.add_argument("--output-tsv-file", dest='output_tsv_file', type=str, help="Output TSV file.")
    arguments = parser.parse_args()
    return arguments


def run():
    args = parse_args()

    df_support = pd.read_csv(args.support_tsv_file, sep="\t")
    df_support = df_support.loc[
        (df_support['num_molecules'] >= args.min_molecules) &
        (df_support['num_cells'] >= args.min_cells),
        :
    ]
    df_support.to_csv(args.output_tsv_file, sep="\t", index=False)

    retained_transcript_ids = set(df_support['transcript_id'].astype(str).unique())

    with pysam.FastaFile(args.fasta_file) as fa_in, open(args.output_fasta_file, 'w') as fa_out:
        for tid in fa_in.references:
            if str(tid) in retained_transcript_ids:
                seq = fa_in.fetch(tid)
                fa_out.write(f">{tid}\n{seq}\n")

    print(f"{len(retained_transcript_ids)} transcripts retained out of {len(df_support)} total.", flush=True)


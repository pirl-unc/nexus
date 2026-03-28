"""
The purpose of this python3 file is to filter RNAbloom2 transcripts.
"""


import argparse
import gzip
import pandas as pd
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description="Filter RNA-Bloom2 assembled transcripts.")
    parser.add_argument("--assembly4-pol-fasta-file", dest='assembly4_pol_fasta_file', type=str, help="RNA-Bloom2 assembly4.pol.fa file.")
    parser.add_argument("--assembly3-map-paf-file", dest='assembly3_map_paf_file', type=str, help="RNA-Bloom2 assembly3.map.paf.gz file.")
    parser.add_argument("--min-mapping-quality", dest='min_mapping_quality', default=30, type=int, help="Minimum mapping quality (default: 30).")
    parser.add_argument("--min-read-support", dest='min_read_support', default=3, type=int, help="Minimum read support (default: 3).")
    parser.add_argument("--min-fraction-match", dest='min_fraction_match', default=0.5, type=float, help="Minimum fraction match (default: 0.5).")
    parser.add_argument("--output-reads-tsv-file", dest='output_reads_tsv_file', type=str, help="Output reads TSV file.")
    parser.add_argument("--output-transcripts-tsv-file", dest='output_transcripts_tsv_file', type=str, help="Output transcripts TSV file.")
    parser.add_argument("--output-fasta-file", dest='output_fasta_file', type=str, help="Output FASTA file.")
    arguments = parser.parse_args()
    return arguments


def run():
    args = parse_args()

    # Step 1. Load the assembled transcript sequences
    # transcripts:
    # {
    #     transcript_id_1: transcript_sequence_1,
    #     ...
    # }
    print("Loading the assembled transcripts FASTA file")
    transcripts = {}
    with open(args.assembly4_pol_fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            transcripts[int(record.id)] = str(record.seq)

    # Step 2. Load the PAF file
    print("Loading the PAF file")
    data = {
        'read_name': [],
        'read_length': [],
        'read_start': [],
        'read_end': [],
        'strand': [],
        'transcript_id': [],
        'transcript_length': [],
        'transcript_start': [],
        'transcript_end': [],
        'num_residue_matches': [],
        'frac_residue_matches': [],
        'alignment_block_length': [],
        'mapping_quality': []
    }
    with gzip.open(args.assembly3_map_paf_file, 'rt') as file:
        for line in file:
            values = line.strip().split('\t')
            frac_residue_matches = float(values[9]) / float(values[6])
            data['read_name'].append(str(values[0]))
            data['read_length'].append(int(values[1]))
            data['read_start'].append(int(values[2]))
            data['read_end'].append(int(values[3]))
            data['strand'].append(str(values[4]))
            data['transcript_id'].append(int(values[5]))
            data['transcript_length'].append(int(values[6]))
            data['transcript_start'].append(int(values[7]))
            data['transcript_end'].append(int(values[8]))
            data['num_residue_matches'].append(int(values[9]))
            data['frac_residue_matches'].append(frac_residue_matches)
            data['alignment_block_length'].append(int(values[10]))
            data['mapping_quality'].append(int(values[11]))
    df_paf = pd.DataFrame(data)
    df_paf_filtered = df_paf[
        (df_paf['mapping_quality'] >= args.min_mapping_quality) &
        (df_paf['frac_residue_matches'] >= args.min_fraction_match)
    ]
    df_grouped = df_paf_filtered.groupby("transcript_id")["read_name"].nunique().reset_index()
    df_grouped = df_grouped.rename(columns={"read_name": "num_read_support"})
    df_grouped = df_grouped[df_grouped['num_read_support'] >= args.min_read_support]
    df_paf_filtered = df_paf_filtered[df_paf_filtered['transcript_id'].isin(df_grouped['transcript_id'].unique())]
    print('%i unique transcripts before filtering' % len(df_paf))
    print('%i unique transcripts after filtering' % len(df_grouped['transcript_id'].unique()))

    # Step 3. Output to TSV files
    df_paf_filtered.to_csv(args.output_reads_tsv_file, index=False, sep='\t')
    data = {
        'transcript_id': [],
        'transcript_sequence': [],
        'num_read_support': []
    }
    for transcript_id in df_grouped['transcript_id'].unique():
        transcript_sequence = transcripts[transcript_id]
        num_read_support = int(
            df_grouped.loc[df_grouped['transcript_id'] == transcript_id, 'num_read_support'].values[0])
        data['transcript_id'].append(transcript_id)
        data['transcript_sequence'].append(transcript_sequence)
        data['num_read_support'].append(num_read_support)
    df_transcripts = pd.DataFrame(data)
    df_transcripts.to_csv(args.output_transcripts_tsv_file, index=False, sep='\t')

    # Step 4. Output to FASTA file
    with open(args.output_fasta_file, 'w') as file:
        for transcript_id in df_grouped['transcript_id'].unique():
            file.write('>%s\n' % transcript_id)
            file.write('%s\n' % transcripts[transcript_id])

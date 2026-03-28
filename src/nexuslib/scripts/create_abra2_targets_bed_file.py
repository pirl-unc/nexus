"""
The purpose of this python3 file is to create an ABRA2 targets BED file.
"""


import pandas as pd
import os
import argparse


def parse_args():
    arg_parser = argparse.ArgumentParser(
        description="""
        Creates an ABRA2 targets BED file based on a GENCODE GTF file.
        The following tasks are performed:
            1. Fetch all exon regions.
            2. Drop duplicates and sort the exon regions.
            3. Output the DataFrame to a BED file.
            4. Merge regions by bedtools.
        """
    )
    arg_parser.add_argument(
        "--gencode_gtf_file",
        dest="gencode_gtf_file",
        type=str,
        required=True,
        help="Input GENCODE GTF file including path (e.g. /<path>/gencode_<version>.gtf)."
    )
    arg_parser.add_argument(
        "--bedtools",
        dest="bedtools",
        type=str,
        required=True,
        help="bedtools path (e.g. /<path>/bedtools.static.binary)."
    )
    arg_parser.add_argument(
        "--output_temp_bed_file",
        dest="output_temp_bed_file",
        type=str,
        required=True,
        help="Output temp BED file including path (e.g. /<path>/gencode_<version>_abra2_targets_temp.bed)."
    )
    arg_parser.add_argument(
        "--output_bed_file",
        dest="output_bed_file",
        type=str,
        required=True,
        help="""
             Output BED file including path (e.g. /<path>/gencode_<version>_abra2_targets.bed). 
             This is the final BED file that should be used as input for ABRA2.
             """
    )
    args = arg_parser.parse_args()
    return args


def run():
    args = parse_args()

    # Step 1. Read the GENCODE GTF file
    df_gencode = pd.read_csv(args.gencode_gtf_file,
                             sep='\t',
                             header=None,
                             comment='#')

    # Step 2. Fetch all exon regions
    data = {
        'Chromosome': [],
        'Start': [],
        'End': []
    }
    for index, row in df_gencode.iterrows():
        curr_chromosome = str(row[0])
        curr_annotation_source = str(row[1])
        curr_type = str(row[2])
        curr_start = int(row[3])
        curr_end = int(row[4])
        if curr_type == 'exon':
            data['Chromosome'].append(curr_chromosome)
            data['Start'].append(curr_start)
            data['End'].append(curr_end)

    # Step 3. Drop duplicates and sort
    df_abra2_targets = pd.DataFrame(data)
    df_abra2_targets.drop_duplicates(inplace=True) # drop duplicates
    df_abra2_targets = df_abra2_targets.groupby(
        ['Chromosome']
    ).apply(
        lambda x: x.sort_values(['Start', 'End'], ascending=True)
    ).reset_index(drop=True)

    # Step 4. Output the dataframe to BED file
    df_abra2_targets.to_csv(args.output_temp_bed_file,
                            header=False,
                            index=False,
                            sep='\t')

    # Step 5. Merge regions by bedtools
    cmd = [
        args.bedtools, 'merge',
        '-i', args.output_temp_bed_file,
        '>', args.output_bed_file
    ]
    cmd = ' '.join(cmd)
    os.system(cmd)

"""
The purpose of this python3 script is to filter unspliced transcripts.
"""


import argparse
import multiprocessing as mp
import pandas as pd
import pysam
from functools import partial
from typing import Dict, List, Set, Tuple
from vstolib.gencode import Gencode


def parse_args():
    parser = argparse.ArgumentParser(description="Filter unspliced transcripts. Transcripts that overlaps any GENCODE transcripts with only 1 exon will be retained.")
    parser.add_argument("--bam-file", dest='bam_file', type=str, help="Input BAM file.")
    parser.add_argument("--fasta-file", dest='fasta_file', type=str, help="Input FASTA file (transcriptome assembly).")
    parser.add_argument("--gencode-gtf-file", dest='gencode_gtf_file', type=str, help="GENCODE GTF file.")
    parser.add_argument("--output-bam-file", dest='output_bam_file', type=str, help="Output BAM file.")
    parser.add_argument("--output-fasta-file", dest='output_fasta_file', type=str, help="Output FASTA file.")
    parser.add_argument("--output-tsv-file", dest='output_tsv_file', type=str, help="Output TSV file.")
    parser.add_argument("--min-mapping-quality", dest='min_mapping_quality', default=30, type=int, help="Minimum mapping quality (default: 30).")
    parser.add_argument("--num-processes", dest='num_processes', default=4, type=int, help="Number of processes (default: 4).")
    arguments = parser.parse_args()
    return arguments


def chunk_dict(dictionary: Dict, num_chunks: int) -> List[Dict]:
    """
    Divide a dictionary into approximately equal chunks.

    Parameters:
        dictionary  : Dictionary to be chunked
        num_chunks  : Number of chunks to create

    Returns:
        List of dictionary chunks
    """
    items = list(dictionary.items())
    chunk_size = len(items) // num_chunks
    remainder = len(items) % num_chunks

    chunks = []
    start = 0

    for i in range(num_chunks):
        # Add one extra item to the first 'remainder' chunks
        current_chunk_size = chunk_size + (1 if i < remainder else 0)
        end = start + current_chunk_size

        # Create a dictionary chunk from the items slice
        chunk_dict = dict(items[start:end])
        chunks.append(chunk_dict)
        start = end

    return chunks


def get_read_ids(bam_file: str) -> Set[str]:
    read_ids = set()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for alignment in bam:
            read_ids.add(alignment.query_name)
    return read_ids


def get_read_ids_passing_mapq(bam_file: str, min_mapping_quality: int) -> Set[str]:
    # Step 1. Get the minimum mapping quality for each read ID
    mapping_qualities: Dict[str, int] = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for alignment in bam:
            if alignment.query_name not in mapping_qualities:
                mapping_qualities[alignment.query_name] = alignment.mapping_quality
            else:
                if alignment.mapping_quality < mapping_qualities[alignment.query_name]:
                    mapping_qualities[alignment.query_name] = alignment.mapping_quality

    # Step 2. Identify reads to filter out
    read_ids_to_keep = set()
    for read_id, mapping_quality in mapping_qualities.items():
        if mapping_quality >= min_mapping_quality:
            read_ids_to_keep.add(read_id)

    return read_ids_to_keep


def index_reads(bam_file: str) -> Dict[str, List[Tuple[str, int, int]]]:
    """
    Index read IDs from a BAM file.

    Parameters:
        bam_file    :   BAM file path.

    Returns:
        Dictionary where keys are read ID and values are a list of tuples (chromosome, start, end).
    """
    reads: Dict[str, List[Tuple[str, int, int]]] = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for alignment in bam:
            if alignment.query_name not in reads:
                reads[alignment.query_name] = []
            reads[alignment.query_name].append((alignment.reference_name, alignment.reference_start, alignment.reference_end))
    return reads


def get_single_exon_transcripts(gencode_gtf_file: str) -> pd.DataFrame:
    gencode = Gencode(gtf_file=gencode_gtf_file, version='', species='')
    exon_counts = gencode.df_exons.groupby('transcript_id')['exon_id'].nunique()
    single_exon_transcript_ids = exon_counts[exon_counts == 1].index
    df_transcripts = gencode.df_transcripts[
        gencode.df_transcripts['transcript_id'].isin(single_exon_transcript_ids)
    ]
    return df_transcripts


def process_alignments(bam_file: str, df_transcripts: pd.DataFrame, reads: Dict[str, List[Tuple[str, int, int]]]) -> Set[str]:
    read_ids_to_keep = set()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read_id, positions in reads.items():
            spliced = False
            overlaps_single_exon_transcript = False
            for (chromosome, start, end) in positions:
                records = bam.fetch(chromosome, start, end)
                for record in records:
                    if record.query_name == read_id:
                        # Check if there is splicing
                        cs_tag = record.get_tag("cs")
                        if '~' in cs_tag:
                            spliced = True

                        # Check if it overlaps a single exon transcript
                        df_matched = df_transcripts[
                            (df_transcripts['chromosome'] == chromosome) &
                            (df_transcripts['start'] <= end) &
                            (df_transcripts['end'] >= start)
                        ]
                        if len(df_matched) > 0:
                            overlaps_single_exon_transcript = True
            if spliced or overlaps_single_exon_transcript:
                read_ids_to_keep.add(read_id)
    return read_ids_to_keep


def filter_fasta_file(fasta_file: str, read_ids_to_keep: Set[str], output_fasta_file: str):
    fasta = pysam.FastaFile(fasta_file)
    with open(output_fasta_file, 'w') as outfile:
        for contig in fasta.references:
            if contig in read_ids_to_keep:
                sequence = fasta.fetch(contig)
                outfile.write(f">{contig}\n{sequence}\n")
    fasta.close()


def write_bam_file(input_bam_file: str, output_bam_file: str, read_ids_to_keep: Set[str]):
    with pysam.AlignmentFile(input_bam_file, "rb") as input_bam:
        with pysam.AlignmentFile(output_bam_file, "wb", template=input_bam) as output_bam:
            for alignment in input_bam:
                if alignment.query_name in read_ids_to_keep:
                    output_bam.write(alignment)
    pysam.index(output_bam_file)


def run():
    # Step 1. Parse arguments
    args = parse_args()

    # Step 2: Get single-exon transcripts
    df_transcripts = get_single_exon_transcripts(args.gencode_gtf_file)
    print("Found %i single-exon transcripts" % len(df_transcripts))

    # Step 3. Count read IDs
    read_ids = get_read_ids(bam_file=args.bam_file)
    print('%i read IDs found in the BAM file.' % len(read_ids))

    # Step 4. Identify read IDs to keep based on mapping quality
    read_ids_passing_mapq = get_read_ids_passing_mapq(bam_file=args.bam_file, min_mapping_quality=args.min_mapping_quality)
    print("%i read IDs meet the minimum mapping quality." % len(read_ids_passing_mapq))

    # Step 5. Index read IDs
    reads = index_reads(bam_file=args.bam_file)
    reads_chunks = chunk_dict(reads, args.num_processes)

    # Step 6. Identify read IDs to keep based on splicing and overlap with single-exon transcripts
    process_func = partial(process_alignments, args.bam_file, df_transcripts)
    with mp.Pool(processes=args.num_processes) as pool:
        results = pool.map(process_func, reads_chunks)
    read_ids_with_splicing = set()
    for result_set in results:
        read_ids_with_splicing.update(result_set)
    print("%i read IDs include splicing or overlap with single-exon transcripts" % len(read_ids_with_splicing))

    # Step 7. Identify eligible read IDs
    eligible_read_ids = read_ids_passing_mapq.intersection(read_ids_with_splicing)
    print("%i read IDs pass the minimum mapping quality and include splicing or overlap with single-exon transcripts" % len(read_ids_with_splicing))

    # Step 8. Filter FASTA file
    filter_fasta_file(fasta_file=args.fasta_file, read_ids_to_keep=eligible_read_ids, output_fasta_file=args.output_fasta_file)

    # Step 9. Write output BAM file
    write_bam_file(args.bam_file, args.output_bam_file, eligible_read_ids)

    # Step 10. Write output TSV file
    data = {
        'read_id': [],
        'status': [],
        'description': []
    }
    for read_id in read_ids:
        if read_id in eligible_read_ids:
            data['read_id'].append(read_id)
            data['status'].append('passed')
            data['description'].append('')
        else:
            description = []
            if read_id not in read_ids_passing_mapq:
                description.append('low_mapping_quality')
            if read_id not in read_ids_with_splicing:
                description.append('unspliced')
            data['read_id'].append(read_id)
            data['status'].append('failed')
            data['description'].append(','.join(description))
    df_results = pd.DataFrame(data)
    df_results.to_csv(args.output_tsv_file, sep='\t', index=False)


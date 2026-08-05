import os
import sys
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *


if __name__ == "__main__":
    # Step 1. Load HLA DNA sequences
    hla_fasta_file = "../../../test/data/fasta/hg38_hla_gen.fa"
    pysam.faidx(hla_fasta_file)
    fasta = pysam.FastaFile(hla_fasta_file)

    sequences = []
    for i, contig in enumerate(fasta.references):
        sequences.append(fasta.fetch(contig, 0, fasta.lengths[i]))

    # Step 2. Create FASTQ files
    create_long_read_fastq_file(
        sequences=sequences,
        output_fastq_file='../../../test/data/fastq/nexus-dna-104-normal_long_read.fastq.gz',
        num_reads=[30] * len(sequences),
        read_length=20000,
        stranded=False
    )
    create_paired_end_fastq_files(
        sequences=sequences,
        output_r1_fastq_file='../../../test/data/fastq/nexus-dna-104-normal_paired-end_read_r1.fastq.gz',
        output_r2_fastq_file = '../../../test/data/fastq/nexus-dna-104-normal_paired-end_read_r2.fastq.gz',
        num_reads=[10000] * len(sequences),
        read_length=151,
        min_insert_size=300,
        max_insert_size=600
    )

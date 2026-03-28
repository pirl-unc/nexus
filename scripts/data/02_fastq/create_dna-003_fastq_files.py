import os
import sys
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("../../../test/data/fasta/GRCh38.p14.genome.chr17.fa.gz")

    # Step 2. Fetch TP53 (chr17:7668421-7687490) sequence
    chromosome = 'chr17'
    start = 7668421
    end = 7687490
    length = end - start + 1
    sequence_normal = str(fasta.fetch(chromosome, start - 1, end))

    # Step 3. Create a somatic deletion at position 7675119-7675148
    reference_position_1 = 7675118
    reference_position_2 = 7675148
    local_position_1 = length - (end - reference_position_1) - 1
    local_position_2 = length - (end - reference_position_2) - 1
    sequence_tumor = sequence_normal
    sequence_tumor = sequence_tumor[:local_position_1 + 1] + sequence_tumor[local_position_2 + 1:]

    # Step 4. Create FASTQ files
    create_long_read_fastq_file(
        sequences=[sequence_tumor, sequence_normal],
        output_fastq_file='../../../test/data/fastq/nexus-dna-003-tumor_long_read.fastq.gz',
        num_reads=[10,10],
        read_length=20000,
        stranded=False
    )
    create_long_read_fastq_file(
        sequences=[sequence_normal],
        output_fastq_file='../../../test/data/fastq/nexus-dna-003-normal_long_read.fastq.gz',
        num_reads=[20],
        read_length=20000,
        stranded=False
    )
    create_paired_end_fastq_files(
        sequences=[sequence_tumor, sequence_normal],
        output_r1_fastq_file='../../../test/data/fastq/nexus-dna-003-tumor_paired-end_read_r1.fastq.gz',
        output_r2_fastq_file = '../../../test/data/fastq/nexus-dna-003-tumor_paired-end_read_r2.fastq.gz',
        num_reads=[5000,5000],
        read_length=151,
        min_insert_size=300,
        max_insert_size=600
    )
    create_paired_end_fastq_files(
        sequences=[sequence_normal],
        output_r1_fastq_file='../../../test/data/fastq/nexus-dna-003-normal_paired-end_read_r1.fastq.gz',
        output_r2_fastq_file = '../../../test/data/fastq/nexus-dna-003-normal_paired-end_read_r2.fastq.gz',
        num_reads=[10000],
        read_length=151,
        min_insert_size=300,
        max_insert_size=600
    )
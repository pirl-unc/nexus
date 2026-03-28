import os
import sys
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("../../../test/data/fasta/GRCh38.p14.genome.chr17.fa.gz")

    # Step 2. Fetch SPEM2 (chr17:7425616-7427568) sequence
    chromosome = 'chr17'
    start = 7425616
    end = 7427568
    length = end - start + 1
    sequence_normal = str(fasta.fetch(chromosome, start - 1, end))

    # Step 3. Create a somatic SNV at position 7426033 (A>G)
    reference_position = 7426033
    local_position = length - (end - reference_position) - 1
    sequence_tumor = sequence_normal
    assert sequence_tumor[local_position] == 'A'
    sequence_tumor = sequence_tumor[:local_position] + 'G' + sequence_tumor[local_position+1:]

    # Step 4. Create FASTQ files
    create_long_read_fastq_file(
        sequences=[sequence_tumor, sequence_normal],
        output_fastq_file='../../../test/data/fastq/nexus-dna-101-tumor_long_read.fastq.gz',
        num_reads=[10,10],
        read_length=20000,
        stranded=False
    )
    create_long_read_fastq_file(
        sequences=[sequence_normal],
        output_fastq_file='../../../test/data/fastq/nexus-dna-101-normal_long_read.fastq.gz',
        num_reads=[20],
        read_length=20000,
        stranded=False
    )
    create_paired_end_fastq_files(
        sequences=[sequence_tumor, sequence_normal],
        output_r1_fastq_file='../../../test/data/fastq/nexus-dna-101-tumor_paired-end_read_r1.fastq.gz',
        output_r2_fastq_file = '../../../test/data/fastq/nexus-dna-101-tumor_paired-end_read_r2.fastq.gz',
        num_reads=[1000,1000],
        read_length=151,
        min_insert_size=300,
        max_insert_size=600
    )
    create_paired_end_fastq_files(
        sequences=[sequence_normal],
        output_r1_fastq_file='../../../test/data/fastq/nexus-dna-101-normal_paired-end_read_r1.fastq.gz',
        output_r2_fastq_file = '../../../test/data/fastq/nexus-dna-101-normal_paired-end_read_r2.fastq.gz',
        num_reads=[2000],
        read_length=151,
        min_insert_size=300,
        max_insert_size=600
    )
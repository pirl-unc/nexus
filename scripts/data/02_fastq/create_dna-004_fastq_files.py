import os
import sys
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("../../../test/data/fasta/GRCh38.p14.genome.chr17.fa.gz")

    # Step 2. Fetch SPNS2 (chr17:4498881-4539035) sequence
    chromosome = 'chr17'
    start_1 = 4498881
    end_1 = 4539035
    length_1 = end_1 - start_1 + 1
    spns2_sequence_normal = str(fasta.fetch(chromosome, start_1 - 1, end_1))

    # Step 3. Fetch PIMREG (chr17:6444455-6451469) sequence
    chromosome = 'chr17'
    start_2 = 6444455
    end_2 = 6451469
    length_2 = end_2 - start_2 + 1
    pimreg_sequence_normal = str(fasta.fetch(chromosome, start_2 - 1, end_2))

    # Step 4. Create a somatic deletion at position 4530700-6447600
    reference_position_1 = 4530700
    reference_position_2 = 6447600
    local_position_1 = length_1 - (end_1 - reference_position_1) - 1
    local_position_2 = length_2 - (end_2 - reference_position_2) - 1
    tumor_sequence = spns2_sequence_normal[:local_position_1] + pimreg_sequence_normal[local_position_2:]

    # Step 5. Create FASTQ files
    create_long_read_fastq_file(
        sequences=[tumor_sequence, spns2_sequence_normal, pimreg_sequence_normal],
        output_fastq_file='../../../test/data/fastq/nexus-dna-004-tumor_long_read.fastq.gz',
        num_reads=[10,10,10],
        read_length=20000,
        stranded=False
    )
    create_long_read_fastq_file(
        sequences=[spns2_sequence_normal, pimreg_sequence_normal],
        output_fastq_file='../../../test/data/fastq/nexus-dna-004-normal_long_read.fastq.gz',
        num_reads=[20, 20],
        read_length=20000,
        stranded=False
    )
    create_paired_end_fastq_files(
        sequences=[tumor_sequence, spns2_sequence_normal, pimreg_sequence_normal],
        output_r1_fastq_file='../../../test/data/fastq/nexus-dna-004-tumor_paired-end_read_r1.fastq.gz',
        output_r2_fastq_file = '../../../test/data/fastq/nexus-dna-004-tumor_paired-end_read_r2.fastq.gz',
        num_reads=[10000,10000,10000],
        read_length=151,
        min_insert_size=300,
        max_insert_size=600
    )
    create_paired_end_fastq_files(
        sequences=[spns2_sequence_normal, pimreg_sequence_normal],
        output_r1_fastq_file='../../../test/data/fastq/nexus-dna-004-normal_paired-end_read_r1.fastq.gz',
        output_r2_fastq_file = '../../../test/data/fastq/nexus-dna-004-normal_paired-end_read_r2.fastq.gz',
        num_reads=[10000, 10000],
        read_length=151,
        min_insert_size=300,
        max_insert_size=600
    )
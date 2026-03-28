import os
import sys
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("../../../test/data/fasta/GRCh38.p14.genome.chr21-22.fa.gz")

    # Step 2. Create tumor and normal DNA
    tumor_dna = []
    normal_dna = []
    for i,chromosome in enumerate(fasta.references):
        length = fasta.lengths[0]
        sequence = fasta.fetch(chromosome, 0, length)

        # Normal DNA
        haplotype_1_sequence_normal = list(str(sequence))
        haplotype_2_sequence_normal = list(str(sequence))
        haplotype_1_snp_positions = range(100000, len(sequence), 1000)
        haplotype_2_snp_positions = range(100000, len(sequence), 1000)
        for snp_position in haplotype_1_snp_positions:
            haplotype_1_sequence_normal[snp_position] = 'A'
        for snp_position in haplotype_2_snp_positions:
            haplotype_2_sequence_normal[snp_position] = 'C'

        # Tumor DNA
        haplotype_1_sequence_tumor = haplotype_1_sequence_normal[0:20000000] + haplotype_1_sequence_normal[30000000:]
        haplotype_2_sequence_tumor = haplotype_2_sequence_normal[0:20000000] + haplotype_2_sequence_normal[30000000:]
        haplotype_1_snv_positions = range(101234, len(haplotype_1_sequence_tumor), 3000)
        haplotype_2_snv_positions = range(105678, len(haplotype_2_sequence_tumor), 3000)
        for snv_position in haplotype_1_snv_positions:
            haplotype_1_sequence_tumor[snv_position] = 'T'
        for snv_position in haplotype_2_snv_positions:
            haplotype_2_sequence_tumor[snv_position] = 'G'

        tumor_dna.append(''.join(haplotype_1_sequence_tumor))
        tumor_dna.append(''.join(haplotype_2_sequence_tumor))
        tumor_dna.append(''.join(haplotype_1_sequence_normal))
        tumor_dna.append(''.join(haplotype_2_sequence_normal))
        normal_dna.append(''.join(haplotype_1_sequence_normal))
        normal_dna.append(''.join(haplotype_2_sequence_normal))

    # Step 3. Create FASTQ files
    create_long_read_fastq_file(
        sequences=tumor_dna,
        output_fastq_file='../../../test/data/fastq/nexus-dna-005-tumor_long_read.fastq.gz',
        num_reads=[1000] * len(tumor_dna),
        read_length=20000,
        stranded=False
    )
    create_long_read_fastq_file(
        sequences=normal_dna,
        output_fastq_file='../../../test/data/fastq/nexus-dna-005-normal_long_read.fastq.gz',
        num_reads=[1000] * len(normal_dna),
        read_length=20000,
        stranded=False
    )
    create_paired_end_fastq_files(
        sequences=tumor_dna,
        output_r1_fastq_file='../../../test/data/fastq/nexus-dna-005-tumor_paired-end_read_r1.fastq.gz',
        output_r2_fastq_file = '../../../test/data/fastq/nexus-dna-005-tumor_paired-end_read_r2.fastq.gz',
        num_reads=[80000] * len(tumor_dna),
        read_length=151,
        min_insert_size=300,
        max_insert_size=600
    )
    create_paired_end_fastq_files(
        sequences=normal_dna,
        output_r1_fastq_file='../../../test/data/fastq/nexus-dna-005-normal_paired-end_read_r1.fastq.gz',
        output_r2_fastq_file = '../../../test/data/fastq/nexus-dna-005-normal_paired-end_read_r2.fastq.gz',
        num_reads=[80000] * len(normal_dna),
        read_length=151,
        min_insert_size=300,
        max_insert_size=600
    )

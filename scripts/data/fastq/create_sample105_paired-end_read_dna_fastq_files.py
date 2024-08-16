import os
import pysam
import random
from typing import List


def reverse_complement(sequence):
    complement_dict = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'a': 't',
        't': 'a',
        'c': 'g',
        'g': 'c',
        'N': 'N',
        'n': 'n'
    }
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(complement_dict[base] for base in reverse_sequence)
    return reverse_complement_sequence

def write_fastq_files(
        sequences: List[str],
        fastq_r1_file: str,
        fastq_r2_file: str,
        num_reads: int,
        read_length: int = 150,
        min_insert_size: int = 300,
        max_insert_size: int = 600
):
    r1_file = open(fastq_r1_file, 'w')
    r2_file = open(fastq_r2_file, 'w')
    for sequence in sequences:
        for i in range(0, num_reads):
            # Read ID
            # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>
            r1_read_id = '@A:1:A:1:1:%i:2 1:N:0:G+A' % (i + 1)
            r2_read_id = '@A:1:A:1:1:%i:2 2:N:0:G+A' % (i + 1)
            while True:
                # Pick a random subsequence
                insert_size = random.randint(min_insert_size, max_insert_size)
                inner_distance = insert_size - (read_length * 2)

                # Pick a random r1 start position
                r1_start = random.randint(0, len(sequence) - insert_size)

                # Read 150bp
                r1_seq = sequence[r1_start:r1_start + read_length]
                r1_seq = ''.join(r1_seq)

                # Determine r2 start position
                r2_start = r1_start + read_length + inner_distance

                # Read 150bp
                r2_seq = sequence[r2_start:r2_start + read_length]
                r2_seq = ''.join(r2_seq)
                r2_seq = reverse_complement(sequence=r2_seq)

                if ('n' in r1_seq) or ('n' in r2_seq) or ('N' in r1_seq) or ('N' in r2_seq):
                    continue

                # Sample base quality scores
                r1_base_quality_scores = [chr(60) for _ in range(0, len(r1_seq))]
                r1_base_quality_scores = ''.join(r1_base_quality_scores)
                r2_base_quality_scores = [chr(60) for _ in range(0, len(r2_seq))]
                r2_base_quality_scores = ''.join(r2_base_quality_scores)

                # Write to R1 file
                r1_file.write(r1_read_id + '\n')
                r1_file.write(r1_seq + '\n')
                r1_file.write('+\n')
                r1_file.write(r1_base_quality_scores + '\n')

                # Write to R2 file
                r2_file.write(r2_read_id + '\n')
                r2_file.write(r2_seq + '\n')
                r2_file.write('+\n')
                r2_file.write(r2_base_quality_scores + '\n')
                break

    r1_file.close()
    r2_file.close()
    os.system('gzip %s' % fastq_r1_file)
    os.system('gzip %s' % fastq_r2_file)


if __name__ == "__main__":
    fasta = pysam.FastaFile("../../../test/data/fasta/hg38_chr21_chr22.fa.gz")
    tumor_dna = []
    normal_dna = []
    for i,chromosome in enumerate(fasta.references):
        sequence = fasta.fetch(chromosome, 0, fasta.lengths[i])

        # Normal DNA
        haplotype_1_sequence_normal = list(str(sequence))
        haplotype_2_sequence_normal = list(str(sequence))
        haplotype_1_snp_positions = range(100000,fasta.lengths[i],1000)
        haplotype_2_snp_positions = range(100000,fasta.lengths[i],1000)
        for snp_position in haplotype_1_snp_positions:
            haplotype_1_sequence_normal[snp_position] = 'A'
        for snp_position in haplotype_2_snp_positions:
            haplotype_2_sequence_normal[snp_position] = 'C'

        # Tumor DNA
        haplotype_1_sequence_tumor = haplotype_1_sequence_normal[0:20000000] + haplotype_1_sequence_normal[30000000:]
        haplotype_2_sequence_tumor = haplotype_2_sequence_normal[0:20000000] + haplotype_2_sequence_normal[30000000:]
        haplotype_1_snv_positions = range(101234,len(haplotype_1_sequence_tumor),3000)
        haplotype_2_snv_positions = range(105678,len(haplotype_2_sequence_tumor),3000)
        for snv_position in haplotype_1_snv_positions:
            haplotype_1_sequence_tumor[snv_position] = 'T'
        for snv_position in haplotype_2_snv_positions:
            haplotype_2_sequence_tumor[snv_position] = 'G'

        tumor_dna.append(''.join(haplotype_1_sequence_tumor))
        tumor_dna.append(''.join(haplotype_2_sequence_tumor))
        normal_dna.append(''.join(haplotype_1_sequence_normal))
        normal_dna.append(''.join(haplotype_2_sequence_normal))

    write_fastq_files(
        sequences=tumor_dna,
        fastq_r1_file='../../../test/data/fastq/sample105tumor_paired-end_read_dna.r1.fastq',
        fastq_r2_file='../../../test/data/fastq/sample105tumor_paired-end_read_dna.r2.fastq',
        num_reads=100000
    )

    write_fastq_files(
        sequences=normal_dna,
        fastq_r1_file='../../../test/data/fastq/sample105normal_paired-end_read_dna.r1.fastq',
        fastq_r2_file='../../../test/data/fastq/sample105normal_paired-end_read_dna.r2.fastq',
        num_reads=100000
    )

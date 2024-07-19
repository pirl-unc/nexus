import os
import pysam
import random


READ_LENGTH = 150
MIN_INSERT_SIZE = 300
MAX_INSERT_SIZE = 600


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


if __name__ == "__main__":
    fasta = pysam.FastaFile("../../../test/data/fasta/hg38_chr21-22.fa.gz")
    r1_file = open('../../../test/data/fastq/hg38_chr21-22_tumor_paired_end_dna.r1.fastq', 'w')
    r2_file = open('../../../test/data/fastq/hg38_chr21-22_tumor_paired_end_dna.r2.fastq', 'w')

    for i,chromosome in enumerate(fasta.references):
        sequence = fasta.fetch(chromosome, 0, fasta.lengths[i])
        haplotype_1_sequence = list(str(sequence))
        haplotype_2_sequence = list(str(sequence))
        haplotype_1_snp_positions = range(5000000,fasta.lengths[i],10000)
        haplotype_2_snp_positions = range(5000000,fasta.lengths[i],10000)
        for snp_position in haplotype_1_snp_positions:
            haplotype_1_sequence[snp_position] = 'A'
        for snp_position in haplotype_2_snp_positions:
            haplotype_2_sequence[snp_position] = 'C'

        # Somatic mutations
        haplotype_1_sequence = haplotype_1_sequence[0:20000000] + haplotype_1_sequence[30000000:]
        haplotype_2_sequence = haplotype_2_sequence[0:20000000] + haplotype_2_sequence[30000000:]
        haplotype_1_snv_positions = range(5001234,len(haplotype_1_sequence),100000)
        haplotype_2_snv_positions = range(5005678,len(haplotype_2_sequence),120000)
        for snv_position in haplotype_1_snv_positions:
            haplotype_1_sequence[snv_position] = 'T'
        for snv_position in haplotype_2_snv_positions:
            haplotype_2_sequence[snv_position] = 'G'

        for i in range(0, 60000):
            # Read ID
            # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>
            r1_read_id = '@A:001:ABC:1:10:%i:20 1:N:0:GA+CA' % (i + 1)
            r2_read_id = '@A:001:ABC:1:10:%i:20 2:N:0:GA+CA' % (i + 1)
            while True:
                if random.randint(0,1) == 0:
                    sequence = haplotype_1_sequence
                else:
                    sequence = haplotype_2_sequence

                # Pick a random subsequence
                insert_size = random.randint(MIN_INSERT_SIZE, MAX_INSERT_SIZE)
                inner_distance = insert_size - (READ_LENGTH * 2)

                # Pick a random r1 start position
                r1_start = random.randint(0, len(sequence) - insert_size)

                # Read 150bp
                r1_seq = sequence[r1_start:r1_start+READ_LENGTH]
                r1_seq = ''.join(r1_seq)

                # Determine r2 start position
                r2_start = r1_start + READ_LENGTH + inner_distance

                # Read 150bp
                r2_seq = sequence[r2_start:r2_start+READ_LENGTH]
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

    os.system('gzip ../../../test/data/fastq/hg38_chr21-22_tumor_paired_end_dna.r1.fastq')
    os.system('gzip ../../../test/data/fastq/hg38_chr21-22_tumor_paired_end_dna.r2.fastq')

import pysam
import random


READ_LENGTH = 150
MIN_INSERT_SIZE = 300
MAX_INSERT_SIZE = 600


def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(complement_dict[base] for base in reverse_sequence)
    return reverse_complement_sequence


if __name__ == "__main__":
    fasta = pysam.FastaFile("data/hg38_tp53_dna.fna")
    sequence = fasta.fetch("NC_000017.11:c7687490-7668421", 0, fasta.lengths[0])
    sequence = list(str(sequence))

    # Point mutations
    sequence[99] = 'A'  # 100 G>A
    sequence[199] = 'A' # 200 G>A
    sequence[299] = 'A' # 300 G>A
    sequence[399] = 'C' # 400 T>C
    sequence[499] = 'T' # 500 A>T

    r1_file = open('hg38_tp53_normal_paired_end_dna.r1.fastq', 'w')
    r2_file = open('hg38_tp53_normal_paired_end_dna.r2.fastq', 'w')

    for i in range(0, 2000):
        # Read ID
        # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>
        r1_read_id = '@A12345:001:ABCDEFGH1:1:1000:%i:2000 1:N:0:GACTGAGT+CACTATCA' % (i + 1)
        r2_read_id = '@A12345:001:ABCDEFGH1:1:1000:%i:2000 2:N:0:GACTGAGT+CACTATCA' % (i + 1)

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

    r1_file.close()
    r2_file.close()
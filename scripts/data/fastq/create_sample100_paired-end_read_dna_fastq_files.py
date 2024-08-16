import os
import pysam
import random


READ_LENGTH = 150
MIN_INSERT_SIZE = 300
MAX_INSERT_SIZE = 600


def create_fastq_files(
        sequence: str,
        output_r1_fastq_file: str,
        output_r2_fastq_file: str,
        num_reads: int
):
    """
    Create FASTQ files.

    Parameters:
        sequence                :   Sequence.
        output_r1_fastq_file    :   Output R1 FASTQ file.
        output_r2_fastq_file    :   Output R2 FASTQ file.
        num_reads               :   Number of reads.
    """
    with open(output_r1_fastq_file, 'w') as f1, open(output_r2_fastq_file, 'w') as f2:
        for i in range(0, num_reads):
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
            f1.write(r1_read_id + '\n')
            f1.write(r1_seq + '\n')
            f1.write('+\n')
            f1.write(r1_base_quality_scores + '\n')

            # Write to R2 file
            f2.write(r2_read_id + '\n')
            f2.write(r2_seq + '\n')
            f2.write('+\n')
            f2.write(r2_base_quality_scores + '\n')

    os.system('gzip %s' % output_r1_fastq_file)
    os.system('gzip %s' % output_r2_fastq_file)


def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(complement_dict[base] for base in reverse_sequence)
    return reverse_complement_sequence


if __name__ == "__main__":
    fasta = pysam.FastaFile("../../../test/data/fasta/hg38_tp53_dna.fa")
    sequence = fasta.fetch("NC_000017.11:c7687490-7668421", 0, fasta.lengths[0])
    sequence = list(str(sequence))

    # Step 1. Tumor
    # Point mutations
    sequence[99] = 'A'  # 100 G>A
    sequence[199] = 'A' # 200 G>A
    sequence[299] = 'A' # 300 G>A
    sequence[399] = 'C' # 400 T>C
    sequence[499] = 'T' # 500 A>T

    # Deletions
    sequence = sequence[0:2000] + sequence[2099:]
    sequence = sequence[0:3000] + sequence[3099:]
    sequence = sequence[0:4000] + sequence[4099:]

    # Insertions
    sequence = sequence[0:10000] + ['ATAGTTACGATTACTGATCGGGGCCCCCTATATATATATCTCTACAAAAAAGCGCTAGCT'] + sequence[10000:]
    sequence = sequence[0:11000] + ['TTTATCTATGTACGTGATCGTAGCTGTAGCTATATGCTAGTCGTAGCTGTAGCTGTACGT'] + sequence[11000:]
    sequence = sequence[0:12000] + ['GAGAGCTAGTGATCTGTAGCTGATTATATGCATCGTAGCTAGCTATCGTATCGTACGTAT'] + sequence[12000:]

    sequence = ''.join(sequence)

    create_fastq_files(
        sequence=sequence,
        output_r1_fastq_file='../../../test/data/fastq/sample100tumor_paired-end_read_dna.r1.fastq',
        output_r2_fastq_file='../../../test/data/fastq/sample100tumor_paired-end_read_dna.r2.fastq',
        num_reads=10000
    )

    # Step 2. Normal
    sequence = fasta.fetch("NC_000017.11:c7687490-7668421", 0, fasta.lengths[0])
    sequence = list(str(sequence))

    # Deletion
    sequence = sequence[0:2000] + sequence[2099:]

    sequence = ''.join(sequence)

    create_fastq_files(
        sequence=sequence,
        output_r1_fastq_file='../../../test/data/fastq/sample100normal_paired-end_read_dna.r1.fastq',
        output_r2_fastq_file='../../../test/data/fastq/sample100normal_paired-end_read_dna.r2.fastq',
        num_reads=10000
    )
import os
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
    fasta = pysam.FastaFile("../../../test/data/fasta/hg38_hla_alleles.fa")
    r1_file = '../../../test/data/fastq/hg38_hla-abc_paired_end_read_rna.r1.fastq'
    r2_file = '../../../test/data/fastq/hg38_hla-abc_paired_end_read_rna.r2.fastq'
    with open(r1_file, 'w') as f1, open(r2_file, 'w') as f2:
        for i,contig in enumerate(fasta.references):
            sequence = fasta.fetch(contig, 0, fasta.lengths[i])
            sequence = list(str(sequence))
            sequence_ = sequence

            for i in range(0, 100):
                # Read ID
                # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>
                r1_read_id = '@A12345:001:ABCDEFGH1:1:1000:%i:2000 1:N:0:GACTGAGT+CACTATCA' % (i + 1)
                r2_read_id = '@A12345:001:ABCDEFGH1:1:1000:%i:2000 2:N:0:GACTGAGT+CACTATCA' % (i + 1)

                # Pick a random subsequence
                insert_size = random.randint(MIN_INSERT_SIZE, MAX_INSERT_SIZE)
                inner_distance = insert_size - (READ_LENGTH * 2)

                # Pick a random r1 start position
                r1_start = random.randint(0, len(sequence_) - insert_size)

                # Read 150bp
                r1_seq = sequence_[r1_start:r1_start+READ_LENGTH]
                r1_seq = ''.join(r1_seq)

                # Determine r2 start position
                r2_start = r1_start + READ_LENGTH + inner_distance

                # Read 150bp
                r2_seq = sequence_[r2_start:r2_start+READ_LENGTH]
                r2_seq = ''.join(r2_seq)
                r2_seq = reverse_complement(sequence=r2_seq)

                # Sample base quality scores
                r1_base_quality_scores = [str(chr(33 + 30)) for _ in range(0, len(r1_seq))]
                r1_base_quality_scores = ''.join(r1_base_quality_scores)
                r2_base_quality_scores = [str(chr(33 + 30)) for _ in range(0, len(r2_seq))]
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

    os.system('gzip %s' % r1_file)
    os.system('gzip %s' % r2_file)
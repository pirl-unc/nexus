import os
import pysam
import random
from typing import List


def create_fastq_file(
        sequences: List[str],
        output_fastq_file: str,
        num_reads: List[int]
):
    """
    Creates a FASTQ file.

    Parameters:
        sequences               :   List of sequences.
        output_fastq_file       :   Output FASTQ file.
        num_reads               :   Number of reads.
    """
    with open(output_fastq_file, 'w') as f:
        for i,sequence in enumerate(sequences):
            for j in range(0, num_reads[i]):
                read_id = '@m64012_%i_%i/%i/ccs' % (random.randint(100000,999999),
                                                    random.randint(100000,999999),
                                                    j+1)
                base_quality_scores = [chr(96) for _ in range(0,len(sequence))]
                base_quality_scores = ''.join(base_quality_scores)
                f.write(read_id + '\n')
                f.write(sequence + '\n')
                f.write('+\n')
                f.write(base_quality_scores + '\n')
    os.system('gzip %s' % output_fastq_file)


if __name__ == "__main__":
    fasta_tp53 = pysam.FastaFile("../../../test/data/fasta/hg38_tp53_dna.fa")
    fasta_ube2g1 = pysam.FastaFile("../../../test/data/fasta/hg38_ube2g1_dna.fa")
    fasta_hpv16 = pysam.FastaFile("../../../test/data/fasta/hpv16_dna.fa")

    # Step 1. Tumor
    sequence_tp53 = fasta_tp53.fetch("NC_000017.11:c7687490-7668421", 0, fasta_tp53.lengths[0])
    sequence_tp53 = list(str(sequence_tp53))
    sequence_ube2g1 = fasta_ube2g1.fetch("NC_000017.11:c4366675-4269259", 0, fasta_ube2g1.lengths[0])
    sequence_ube2g1 = list(str(sequence_ube2g1))
    sequence_hpv16 = fasta_hpv16.fetch("K02718.1", 0, fasta_hpv16.lengths[0])
    sequence_hpv16 = list(str(sequence_hpv16))

    # Point mutations
    for i in range(0,len(sequence_tp53),1000):
        sequence_tp53[i] = 'A'
    for i in range(0,len(sequence_ube2g1),1000):
        sequence_ube2g1[i] = 'A'

    # Viral integration
    sequence_translocation = sequence_tp53[0:10000] + list('AACGATCGTAGCTAGTCTAGCTGAAGCGAT') + sequence_hpv16[2000:]
    sequence_translocation = ''.join(sequence_translocation)
    sequence_tp53 = ''.join(sequence_tp53)
    sequence_ube2g1 = ''.join(sequence_ube2g1)
    sequence_hpv16 = ''.join(sequence_hpv16)

    create_fastq_file(
        sequences=[sequence_translocation, sequence_tp53, sequence_ube2g1, sequence_hpv16],
        output_fastq_file='../../../test/data/fastq/sample003tumor_long_read_dna.fastq',
        num_reads=[60,60,60,30]
    )

    # Step 2. Normal
    create_fastq_file(
        sequences=[sequence_tp53, sequence_ube2g1],
        output_fastq_file='../../../test/data/fastq/sample003normal_long_read_dna.fastq',
        num_reads=[60,60]
    )

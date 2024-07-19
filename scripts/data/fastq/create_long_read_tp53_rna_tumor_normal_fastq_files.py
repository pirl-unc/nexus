import os
import pysam
import random
from typing import List


def create_fastq_file(
        sequences: List[str],
        output_fastq_file: str,
        num_reads: int
):
    """
    Creates a FASTQ file.

    Parameters:
        sequences               :   List of sequences.
        output_fastq_file       :   Output FASTQ file.
    """
    with open(output_fastq_file, 'w') as f:
        i = 1
        for sequence in sequences:
            for _ in range(0, num_reads):
                read_id = '@m64012_%i_%i/%i/ccs' % (random.randint(100000,999999),
                                                    random.randint(100000,999999),
                                                    i)
                base_quality_scores = [chr(96) for _ in range(0,len(sequence))]
                base_quality_scores = ''.join(base_quality_scores)
                f.write(read_id + '\n')
                f.write(sequence + '\n')
                f.write('+\n')
                f.write(base_quality_scores + '\n')
                i += 1
    os.system('gzip %s' % output_fastq_file)


if __name__ == "__main__":
    fasta = pysam.FastaFile("../../../test/data/fasta/hg38_tp53_rna.fa")
    ref_sequence = fasta.fetch("ENST00000269305.9", 0, fasta.lengths[0])
    var_sequence = list(str(ref_sequence))

    # Step 1. Tumor
    # Point mutations
    var_sequence[99] = 'A'  # 100 G>A
    var_sequence[199] = 'A' # 200 G>A
    var_sequence[299] = 'A' # 300 G>A
    var_sequence[399] = 'C' # 400 T>C
    var_sequence[499] = 'T' # 500 A>T

    # Deletions
    var_sequence = var_sequence[0:500] + var_sequence[599:]
    var_sequence = var_sequence[0:800] + var_sequence[899:]
    var_sequence = var_sequence[0:1100] + var_sequence[1199:]

    # Insertions
    var_sequence = var_sequence[0:1500] + ['ATAGTTACGATTACTGATCGGGGCCCCCTATATATATATCTCTACAAAAAAGCGCTAGCT'] + var_sequence[1500:]
    var_sequence = var_sequence[0:1800] + ['TTTATCTATGTACGTGATCGTAGCTGTAGCTATATGCTAGTCGTAGCTGTAGCTGTACGT'] + var_sequence[1800:]
    var_sequence = var_sequence[0:2100] + ['GAGAGCTAGTGATCTGTAGCTGATTATATGCATCGTAGCTAGCTATCGTATCGTACGTAT'] + var_sequence[2100:]

    var_sequence = ''.join(var_sequence)
    ref_sequence = ''.join(ref_sequence)

    create_fastq_file(
        sequences=[ref_sequence, var_sequence],
        output_fastq_file='../../../test/data/fastq/hg38_tp53_tumor_long_read_rna.fastq',
        num_reads=200
    )

    # Step 2. Normal
    create_fastq_file(
        sequences=[ref_sequence,],
        output_fastq_file='../../../test/data/fastq/hg38_tp53_normal_long_read_rna.fastq',
        num_reads=200
    )

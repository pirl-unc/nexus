
import pysam
import random


if __name__ == "__main__":
    fasta = pysam.FastaFile("hg38_tp53_ ENST00000269305_rna.fasta")
    sequence = fasta.fetch("NM_000546.6", 0, fasta.lengths[0])
    sequence = list(str(sequence))

    # Point mutations
    sequence[99] = 'A'  # 100 G>A
    sequence[199] = 'A' # 200 G>A
    sequence[299] = 'A' # 300 G>A
    sequence[399] = 'C' # 400 T>C
    sequence[499] = 'T' # 500 A>T

    # Deletions
    sequence = sequence[0:500] + sequence[599:]
    sequence = sequence[0:800] + sequence[899:]
    sequence = sequence[0:1100] + sequence[1199:]

    # Insertions
    sequence = sequence[0:1500] + ['ATAGTTACGATTACTGATCGGGGCCCCCTATATATATATCTCTACAAAAAAGCGCTAGCT'] + sequence[1500:]
    sequence = sequence[0:1800] + ['TTTATCTATGTACGTGATCGTAGCTGTAGCTATATGCTAGTCGTAGCTGTAGCTGTACGT'] + sequence[1800:]
    sequence = sequence[0:2100] + ['GAGAGCTAGTGATCTGTAGCTGATTATATGCATCGTAGCTAGCTATCGTATCGTACGTAT'] + sequence[2100:]

    sequence = ''.join(sequence)

    with open('hg38_tp53_variants_rna.fastq', 'w') as f:
        for i in range(0, 100):
            read_id = '@m64013_%i_%i/%i/ccs' % (random.randint(100000,999999),
                                                random.randint(100000,999999),
                                                i+1)
            base_quality_scores = [chr(96) for _ in range(0,len(sequence))]
            base_quality_scores = ''.join(base_quality_scores)
            f.write(read_id + '\n')
            f.write(sequence + '\n')
            f.write('+\n')
            f.write(base_quality_scores + '\n')

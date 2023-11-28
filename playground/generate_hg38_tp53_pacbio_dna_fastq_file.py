import pysam
import random


if __name__ == "__main__":
    fasta = pysam.FastaFile("hg38_tp53_dna.fna")
    sequence = fasta.fetch("NC_000017.11:c7687490-7668421", 0, fasta.lengths[0])
    sequence = list(str(sequence))

    # Point mutations
    # sequence[99] = 'A'  # 100 G>A
    # sequence[199] = 'A' # 200 G>A
    # sequence[299] = 'A' # 300 G>A
    # sequence[399] = 'C' # 400 T>C
    # sequence[499] = 'T' # 500 A>T
    #
    # # Deletions
    # sequence = sequence[0:2000] + sequence[2099:]
    # sequence = sequence[0:3000] + sequence[3099:]
    # sequence = sequence[0:4000] + sequence[4099:]
    #
    # # Insertions
    # sequence = sequence[0:10000] + ['ATAGTTACGATTACTGATCGGGGCCCCCTATATATATATCTCTACAAAAAAGCGCTAGCT'] + sequence[10000:]
    # sequence = sequence[0:11000] + ['TTTATCTATGTACGTGATCGTAGCTGTAGCTATATGCTAGTCGTAGCTGTAGCTGTACGT'] + sequence[11000:]
    # sequence = sequence[0:12000] + ['GAGAGCTAGTGATCTGTAGCTGATTATATGCATCGTAGCTAGCTATCGTATCGTACGTAT'] + sequence[12000:]

    sequence = ''.join(sequence)

    with open('hg38_tp53_normal_dna.fastq', 'w') as f:
        for i in range(0, 100):
            read_id = '@m64012_%i_%i/%i/ccs' % (random.randint(100000,999999),
                                                random.randint(100000,999999),
                                                i+1)
            base_quality_scores = [chr(96) for _ in range(0,len(sequence))]
            base_quality_scores = ''.join(base_quality_scores)
            f.write(read_id + '\n')
            f.write(sequence + '\n')
            f.write('+\n')
            f.write(base_quality_scores + '\n')

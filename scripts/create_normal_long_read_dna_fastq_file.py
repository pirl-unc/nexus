import pysam
import random


if __name__ == "__main__":
    fasta = pysam.FastaFile("data/hg38_tp53_dna.fna")
    sequence = fasta.fetch("NC_000017.11:c7687490-7668421", 0, fasta.lengths[0])
    sequence = list(str(sequence))
    sequence = ''.join(sequence)

    with open('hg38_tp53_normal_long_read_dna.fastq', 'w') as f:
        for i in range(0, 60):
            read_id = '@m64012_%i_%i/%i/ccs' % (random.randint(100000,999999),
                                                random.randint(100000,999999),
                                                i+1)
            base_quality_scores = [chr(96) for _ in range(0,len(sequence))]
            base_quality_scores = ''.join(base_quality_scores)
            f.write(read_id + '\n')
            f.write(sequence + '\n')
            f.write('+\n')
            f.write(base_quality_scores + '\n')

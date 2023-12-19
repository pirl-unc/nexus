
import pysam
import random


if __name__ == "__main__":
    fasta = pysam.FastaFile("data/hg38_tp53_ ENST00000269305_rna.fasta")
    ref_sequence = fasta.fetch("NM_000546.6", 0, fasta.lengths[0])
    ref_sequence = str(ref_sequence)

    with open('hg38_tp53_normal_long_read_rna.fastq', 'w') as f:
        for i in range(0, 100):
            read_id = '@m64013_%i_%i/%i/ccs' % (random.randint(100000,999999),
                                                random.randint(100000,999999),
                                                i+1)
            base_quality_scores = [chr(96) for _ in range(0,len(ref_sequence))]
            base_quality_scores = ''.join(base_quality_scores)
            f.write(read_id + '\n')
            f.write(ref_sequence + '\n')
            f.write('+\n')
            f.write(base_quality_scores + '\n')

import os
import sys
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *
from vstolib.gencode import Gencode


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("../../../test/data/fasta/GRCh38.p14.genome.chr17.fa.gz")

    # Step 2. Load GENCODE
    gencode = Gencode(
        gtf_file="../../../test/data/gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz",
        version='v41',
        species='human',
        levels=[1,2],
        types=['protein_coding']
    )

    # Step 3. Fetch FGF11 sequence
    df_transcript_fgf11 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000293829']
    df_exons_fgf11 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_fgf11['transcript_id'].values[0]]
    df_exons_fgf11.sort_values(by=['number'], ascending=True, inplace=True) # FGF11 is on the forward strand
    fgf11_tumor_sequence = ''
    fgf11_normal_sequence = ''
    for _,row in df_exons_fgf11.iterrows():
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        if exon_number >= 2 and exon_number <= 4:
            fgf11_tumor_sequence = fgf11_tumor_sequence + normal_sequence
        fgf11_normal_sequence = fgf11_normal_sequence + normal_sequence
    tumor_sequence = fgf11_tumor_sequence + fgf11_tumor_sequence

    # Step 4. Create FASTQ files
    create_long_read_fastq_file(
        sequences=[tumor_sequence, fgf11_normal_sequence],
        output_fastq_file='../../../test/data/fastq/nexus-rna-102-tumor_long_read.fastq.gz',
        num_reads=[10,10],
        read_length=20000,
        stranded=True
    )
    create_paired_end_fastq_files(
        sequences=[tumor_sequence, fgf11_normal_sequence],
        output_r1_fastq_file='../../../test/data/fastq/nexus-rna-102-tumor_paired-end_read_r1.fastq.gz',
        output_r2_fastq_file = '../../../test/data/fastq/nexus-rna-102-tumor_paired-end_read_r2.fastq.gz',
        num_reads=[1000,1000],
        read_length=151,
        min_insert_size=300,
        max_insert_size=600
    )

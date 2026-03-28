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

    # Step 3. Fetch SPNS2 sequence
    df_transcript_spns2 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000329078']
    df_exons_spns2 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_spns2['transcript_id'].values[0]]
    df_exons_spns2.sort_values(by=['number'], ascending=True, inplace=True) # SPNS2 is on the forward strand
    spns2_tumor_sequence = ''
    spns2_normal_sequence = ''
    for _,row in df_exons_spns2.iterrows():
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        if exon_number == 4:
            # Create a DEL 4530700-6447600
            reference_position = 4530700
            local_position = len(normal_sequence) - (end - reference_position) - 1
            tumor_sequence = normal_sequence[:local_position + 1]
        elif exon_number < 4:
            tumor_sequence = normal_sequence
        else:
            tumor_sequence = ''
        spns2_tumor_sequence = spns2_tumor_sequence + tumor_sequence
        spns2_normal_sequence = spns2_normal_sequence + normal_sequence

    # Step 4. Fetch PIMREG sequence
    df_transcript_pimreg = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000250056']
    df_exons_pimreg = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_pimreg['transcript_id'].values[0]]
    df_exons_pimreg.sort_values(by=['number'], ascending=True, inplace=True) # PIMREG is on the forward strand
    pimreg_tumor_sequence = ''
    pimreg_normal_sequence = ''
    for _,row in df_exons_pimreg.iterrows():
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        if exon_number == 3:
            # Create a DEL 4530700-6447600
            reference_position = 6447600
            local_position = len(normal_sequence) - (end - reference_position) - 1
            tumor_sequence = normal_sequence[local_position:]
        elif exon_number > 3:
            tumor_sequence = normal_sequence
        else:
            tumor_sequence = ''
        pimreg_tumor_sequence = pimreg_tumor_sequence + tumor_sequence
        pimreg_normal_sequence = pimreg_normal_sequence + normal_sequence

    # Step 5. Create tumor sequences
    tumor_sequence = spns2_tumor_sequence + pimreg_tumor_sequence

    # Step 6. Create FASTQ files
    create_long_read_fastq_file(
        sequences=[tumor_sequence, spns2_normal_sequence, pimreg_normal_sequence],
        output_fastq_file='../../../test/data/fastq/nexus-rna-004-tumor_long_read.fastq.gz',
        num_reads=[30,30,30],
        read_length=20000,
        stranded=True
    )
    create_paired_end_fastq_files(
        sequences=[tumor_sequence, spns2_normal_sequence, pimreg_normal_sequence],
        output_r1_fastq_file='../../../test/data/fastq/nexus-rna-004-tumor_paired-end_read_r1.fastq.gz',
        output_r2_fastq_file = '../../../test/data/fastq/nexus-rna-004-tumor_paired-end_read_r2.fastq.gz',
        num_reads=[2000,2000,2000],
        read_length=151,
        min_insert_size=300,
        max_insert_size=600
    )

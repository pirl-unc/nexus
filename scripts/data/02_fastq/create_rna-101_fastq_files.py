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

    # Step 3. Fetch SPEM2 sequence and create a somatic SNV at positions:
    # 7425743 (A>G)
    # 7426033 (A>G)
    # 7426653 (A>G)
    df_transcript_spem2 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000574034']
    df_exons_spem2 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_spem2['transcript_id'].values[0]]
    df_exons_spem2.sort_values(by=['number'], ascending=True, inplace=True) # SPEM2 is on the forward strand
    spem2_tumor_sequence = ''
    spem2_normal_sequence = ''
    for _,row in df_exons_spem2.iterrows():
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        if exon_number == 1:
            # Create a SNV at position 7425743
            reference_position = 7425743
            local_position = len(normal_sequence) - (end - reference_position) - 1
            assert normal_sequence[local_position] == 'A'
            tumor_sequence = normal_sequence[:local_position] + 'G' + normal_sequence[local_position+1:]
        elif exon_number == 2:
            # Create a SNV at position 7426033
            reference_position = 7426033
            local_position = len(normal_sequence) - (end - reference_position) - 1
            assert normal_sequence[local_position] == 'A'
            tumor_sequence = normal_sequence[:local_position] + 'G' + normal_sequence[local_position+1:]
        elif exon_number == 3:
            # Create a SNV at position 7426653
            reference_position = 7426653
            local_position = len(normal_sequence) - (end - reference_position) - 1
            assert normal_sequence[local_position] == 'A'
            tumor_sequence = normal_sequence[:local_position] + 'G' + normal_sequence[local_position+1:]
        else:
            tumor_sequence = normal_sequence
        spem2_tumor_sequence = spem2_tumor_sequence + tumor_sequence
        spem2_normal_sequence = spem2_normal_sequence + normal_sequence

    # Step 4. Create FASTQ files
    create_long_read_fastq_file(
        sequences=[spem2_tumor_sequence, spem2_normal_sequence],
        output_fastq_file='../../../test/data/fastq/nexus-rna-101-tumor_long_read.fastq.gz',
        num_reads=[20,20],
        read_length=20000,
        stranded=True
    )
    create_paired_end_fastq_files(
        sequences=[spem2_tumor_sequence, spem2_normal_sequence],
        output_r1_fastq_file='../../../test/data/fastq/nexus-rna-101-tumor_paired-end_read_r1.fastq.gz',
        output_r2_fastq_file = '../../../test/data/fastq/nexus-rna-101-tumor_paired-end_read_r2.fastq.gz',
        num_reads=[2000,4000],
        read_length=151,
        min_insert_size=300,
        max_insert_size=600
    )

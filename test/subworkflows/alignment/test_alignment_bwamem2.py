import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_alignment_bwamem2_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    paired_end_read_tumor_dna_r1_fastq_file = get_data_path(name='fastq/nexus-dna-001-tumor_paired-end_read_r1.fastq.gz')
    paired_end_read_tumor_dna_r2_fastq_file = get_data_path(name='fastq/nexus-dna-001-tumor_paired-end_read_r2.fastq.gz')
    paired_end_read_normal_dna_r1_fastq_file = get_data_path(name='fastq/nexus-dna-001-normal_paired-end_read_r1.fastq.gz')
    paired_end_read_normal_dna_r2_fastq_file = get_data_path(name='fastq/nexus-dna-001-normal_paired-end_read_r2.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    abra2_targets_bed_file = get_data_path(name='indices/abra2/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.abra2_targets.bed')
    known_sites_vcf_file = get_data_path(name='vcf/known_sites_hg38_chr17.vcf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_alignment_bwamem2'
    work_dir = temp_dir + '/work/test_alignment_bwamem2'
    output_dir = temp_dir + '/outputs/test_alignment_bwamem2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor', 'nexus-dna-001-normal'],
        'fastq_file_1': [paired_end_read_tumor_dna_r1_fastq_file, paired_end_read_normal_dna_r1_fastq_file],
        'fastq_file_2': [paired_end_read_tumor_dna_r2_fastq_file, paired_end_read_normal_dna_r2_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--abra2_targets_bed_file', abra2_targets_bed_file,
        '--known_sites_vcf_files', known_sites_vcf_file,
        '--chromosomes', '"chr17"',
        '--perform_local_indel_realignment', 'false',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='alignment_bwamem2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

def test_alignment_bwamem2_2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    paired_end_read_tumor_dna_r1_fastq_file = get_data_path(name='fastq/nexus-dna-001-tumor_paired-end_read_r1.fastq.gz')
    paired_end_read_tumor_dna_r2_fastq_file = get_data_path(name='fastq/nexus-dna-001-tumor_paired-end_read_r2.fastq.gz')
    paired_end_read_normal_dna_r1_fastq_file = get_data_path(name='fastq/nexus-dna-001-normal_paired-end_read_r1.fastq.gz')
    paired_end_read_normal_dna_r2_fastq_file = get_data_path(name='fastq/nexus-dna-001-normal_paired-end_read_r2.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    abra2_targets_bed_file = get_data_path(name='indices/abra2/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.abra2_targets.bed')
    known_sites_vcf_file = get_data_path(name='vcf/known_sites_hg38_chr17.vcf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_alignment_bwamem2_abra2'
    work_dir = temp_dir + '/work/test_alignment_bwamem2_abra2'
    output_dir = temp_dir + '/outputs/test_alignment_bwamem2_abra2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor', 'nexus-dna-001-normal'],
        'fastq_file_1': [paired_end_read_tumor_dna_r1_fastq_file, paired_end_read_normal_dna_r1_fastq_file],
        'fastq_file_2': [paired_end_read_tumor_dna_r2_fastq_file, paired_end_read_normal_dna_r2_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--abra2_targets_bed_file', abra2_targets_bed_file,
        '--known_sites_vcf_files', known_sites_vcf_file,
        '--abra2_temp_dir', os.getcwd() + '/tmp/abra2_temp/',
        '--chromosomes', '"chr17"',
        '--perform_local_indel_realignment', 'true',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='alignment_bwamem2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

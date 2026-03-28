import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_calling_mutect2_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    tumor_bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    normal_bam_file = get_data_path(name='bam/nexus-dna-001-normal-paired-end-read_fixmate_markeddup_recalibrated.bam')
    normal_bam_bai_file = get_data_path(name='bam/nexus-dna-001-normal-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    small_exac_commmon_3_hg38_vcf_file = get_data_path(name='vcf/small_exac_common_3.hg38.vcf')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_mutect2_1'
    work_dir = temp_dir + '/work/test_variant_calling_mutect2_1'
    output_dir = temp_dir + '/outputs/test_variant_calling_mutect2_1'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'tumor_bam_file': [tumor_bam_file],
        'tumor_bam_bai_file': [tumor_bam_bai_file],
        'normal_bam_file': [normal_bam_file],
        'normal_bam_bai_file': [normal_bam_bai_file],
        'normal_sample_id': ['nexus-dna-001-normal-paired-end-read'],
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--mutect2_germline_resource_vcf_file', '""',
        '--mutect2_panel_of_normals_vcf_file', '""',
        '--getpileupsummaries_variant_vcf_file', small_exac_commmon_3_hg38_vcf_file,
        '--chromosomes', '"chr17"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_mutect2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_variant_calling_mutect2_2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    tumor_bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    normal_bam_file = get_data_path(name='bam/nexus-dna-001-normal-paired-end-read_fixmate_markeddup_recalibrated.bam')
    normal_bam_bai_file = get_data_path(name='bam/nexus-dna-001-normal-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_mutect2_2'
    work_dir = temp_dir + '/work/test_variant_calling_mutect2_2'
    output_dir = temp_dir + '/outputs/test_variant_calling_mutect2_2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'tumor_bam_file': [tumor_bam_file],
        'tumor_bam_bai_file': [tumor_bam_bai_file],
        'normal_bam_file': [normal_bam_file],
        'normal_bam_bai_file': [normal_bam_bai_file],
        'normal_sample_id': ['nexus-dna-001-normal-paired-end-read'],
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--mutect2_germline_resource_vcf_file', '""',
        '--mutect2_panel_of_normals_vcf_file', '""',
        '--getpileupsummaries_variant_vcf_file', '""',
        '--chromosomes', '"chr17"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_mutect2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

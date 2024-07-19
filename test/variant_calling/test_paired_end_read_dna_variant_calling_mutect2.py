import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path, get_alias_path


def test_paired_end_read_dna_variant_calling_mutect2_human():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/hg38_tp53_tumor_paired-end_read_dna.bam')
    tumor_bam_bai_file = get_data_path(name='bam/hg38_tp53_tumor_paired-end_read_dna.bam.bai')
    normal_bam_file = get_data_path(name='bam/hg38_tp53_normal_paired-end_read_dna.bam')
    normal_bam_bai_file = get_data_path(name='bam/hg38_tp53_normal_paired-end_read_dna.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa')
    reference_genome_fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa.fai')
    reference_genome_fasta_dict_file = get_data_path(name='fasta/hg38_chr17_1-8000000.dict')
    small_exac_commmon_3_hg38_vcf_file = get_data_path(name='vcf/small_exac_common_3.hg38.vcf')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_paired_end_read_dna_variant_calling_mutect2_human'
    work_dir = temp_dir + '/work/test_paired_end_read_dna_variant_calling_mutect2_human'
    output_dir = temp_dir + '/outputs/test_paired_end_read_dna_variant_calling_mutect2_human'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample001'],
        'tumor_bam_file': [tumor_bam_file],
        'tumor_bam_bai_file': [tumor_bam_bai_file],
        'normal_bam_file': [normal_bam_file],
        'normal_bam_bai_file': [normal_bam_bai_file],
        'tumor_sample_id': ['sample001'],
        'normal_sample_id': ['sample001']
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--reference_genome_fasta_dict_file', reference_genome_fasta_dict_file,
        '--chromosomes', '"chr17"',
        '--is_human', '1',
        '--params_gatk4mutect2', '""',
        '--params_gatk4getpileupsummaries', '"-V %s -L %s"' % (small_exac_commmon_3_hg38_vcf_file, small_exac_commmon_3_hg38_vcf_file),
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='paired-end_read_dna_variant_calling_mutect2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_paired_end_read_dna_variant_calling_mutect2_nonhuman():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/hg38_tp53_tumor_paired-end_read_dna.bam')
    tumor_bam_bai_file = get_data_path(name='bam/hg38_tp53_tumor_paired-end_read_dna.bam.bai')
    normal_bam_file = get_data_path(name='bam/hg38_tp53_normal_paired-end_read_dna.bam')
    normal_bam_bai_file = get_data_path(name='bam/hg38_tp53_normal_paired-end_read_dna.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa')
    reference_genome_fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa.fai')
    reference_genome_fasta_dict_file = get_data_path(name='fasta/hg38_chr17_1-8000000.dict')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_paired_end_read_dna_variant_calling_mutect2_nonhuman'
    work_dir = temp_dir + '/work/test_paired_end_read_dna_variant_calling_mutect2_nonhuman'
    output_dir = temp_dir + '/outputs/test_paired_end_read_dna_variant_calling_mutect2_nonhuman'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample001'],
        'tumor_bam_file': [tumor_bam_file],
        'tumor_bam_bai_file': [tumor_bam_bai_file],
        'normal_bam_file': [normal_bam_file],
        'normal_bam_bai_file': [normal_bam_bai_file],
        'tumor_sample_id': ['sample001'],
        'normal_sample_id': ['sample001']
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--reference_genome_fasta_dict_file', reference_genome_fasta_dict_file,
        '--chromosomes', '"chr17"',
        '--is_human', '0',
        '--params_gatk4mutect2', '""',
        '--params_gatk4getpileupsummaries', '""',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='paired-end_read_dna_variant_calling_mutect2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


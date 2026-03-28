import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_calling_dysgu_somatic_lr():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file_1 = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file_1 = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    bam_file_2 = get_data_path(name='bam/nexus-dna-003-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file_2 = get_data_path(name='bam/nexus-dna-003-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_dysgu_somatic_lr'
    work_dir = temp_dir + '/work/test_variant_calling_dysgu_somatic_lr'
    output_dir = temp_dir + '/outputs/test_variant_calling_dysgu_somatic_lr'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-002-tumor'],
        'tumor_bam_file': [bam_file_1],
        'tumor_bam_bai_file': [bam_bai_file_1],
        'normal_bam_file': [bam_file_2],
        'normal_bam_bai_file': [bam_bai_file_2]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--params_dysgu_run', '"--mode pacbio-revio --min-support 3 --min-size 30 --mq 20"',
        '--params_dysgu_filter', '"--support-fraction 0.05 --min-mapq 20"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_dysgu-somatic.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_variant_calling_dysgu_somatic_pe():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/nexus-dna-002-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    tumor_bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    normal_bam_file = get_data_path(name='bam/nexus-dna-003-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    normal_bam_bai_file = get_data_path(name='bam/nexus-dna-003-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_dysgu_somatic_pe'
    work_dir = temp_dir + '/work/test_variant_calling_dysgu_somatic_pe'
    output_dir = temp_dir + '/outputs/test_variant_calling_dysgu_somatic_pe'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-002-tumor'],
        'tumor_bam_file': [tumor_bam_file],
        'tumor_bam_bai_file': [tumor_bam_bai_file],
        'normal_bam_file': [normal_bam_file],
        'normal_bam_bai_file': [normal_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--params_dysgu_run', '"--mode pe --min-support 3 --min-size 10 --mq 20"',
        '--params_dysgu_filter', '"--support-fraction 0.05 --min-mapq 20"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_dysgu-somatic.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


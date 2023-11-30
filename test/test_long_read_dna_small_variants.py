import pandas as pd
import os
from .data import get_data_path
from nexuslib.main import run_workflow


def test_long_read_dna_small_variants_singularity():
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    long_read_tumor_dna_bam_file = get_data_path(name='hg38_tp53_tumor_long_read_dna_minimap2_mdtagged_sorted.bam')
    long_read_tumor_dna_bam_bai_file = get_data_path(name='hg38_tp53_tumor_long_read_dna_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='hg38_chr17_1-8000000.fa')
    temp_dir = os.getcwd() + '/tmp/'
    intermediate_dir = temp_dir + 'intermediate/test_long_read_dna_small_variants/'
    work_dir = temp_dir + 'work/test_long_read_dna_small_variants/'
    output_dir = temp_dir + 'outputs/test_long_read_dna_small_variants/'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor'],
        'bam_file': [long_read_tumor_dna_bam_file],
        'bam_bai_file': [long_read_tumor_dna_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--containerization', 'singularity',
        '--deepvariant_bin_path', '/opt/deepvariant/bin/run_deepvariant',
        '--deepvariant_input_path', '/home/runner/work/nexus/nexus/',
        '--deepvariant_output_path', '/home/runner/work/nexus/nexus/',
        '--deepvariant_bin_version', '1.6.0',
        '--deepvariant_model_type', 'WGS',
        '--output_dir', output_dir,
        '-w', work_dir
    ]
    run_workflow(workflow='long_read_dna_small_variants.nf',
                 workflow_args=workflow_args)


def test_long_read_dna_small_variants_docker():
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    long_read_tumor_dna_bam_file = get_data_path(name='hg38_tp53_tumor_long_read_dna_minimap2_mdtagged_sorted.bam')
    long_read_tumor_dna_bam_bai_file = get_data_path(name='hg38_tp53_tumor_long_read_dna_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='hg38_chr17_1-8000000.fa')
    temp_dir = os.getcwd() + '/tmp/'
    intermediate_dir = temp_dir + 'intermediate/test_long_read_dna_small_variants/'
    work_dir = temp_dir + 'work/test_long_read_dna_small_variants/'
    output_dir = temp_dir + 'outputs/test_long_read_dna_small_variants/'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor'],
        'bam_file': [long_read_tumor_dna_bam_file],
        'bam_bai_file': [long_read_tumor_dna_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--containerization', 'docker',
        '--deepvariant_bin_path', '/opt/deepvariant/bin/run_deepvariant',
        '--deepvariant_input_path', '/Users/leework/Documents/Research/projects/nexus/',
        '--deepvariant_output_path', '/var/folders/',
        '--deepvariant_bin_version', '1.6.0',
        '--deepvariant_model_type', 'WGS',
        '--output_dir', output_dir,
        '-w', work_dir
    ]
    run_workflow(workflow='long_read_dna_small_variants.nf',
                 workflow_args=workflow_args)
import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_long_read_dna_variant_calling_deepvariant_github():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_dna_bam_file = get_data_path(name='bam/sample001tumor_minimap2_mdtagged_sorted.bam')
    long_read_tumor_dna_bam_bai_file = get_data_path(name='bam/sample001tumor_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_dna_variant_calling_deepvariant'
    work_dir = temp_dir + '/work/test_long_read_dna_variant_calling_deepvariant'
    output_dir = temp_dir + '/outputs/test_long_read_dna_variant_calling_deepvariant'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample001tumor'],
        'bam_file': [long_read_tumor_dna_bam_file],
        'bam_bai_file': [long_read_tumor_dna_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--containerization', 'docker',
        '--deepvariant_input_path', '/home/runner/work/nexus/nexus/',
        '--deepvariant_output_path', '/tmp/',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='long_read_dna_variant_calling_deepvariant.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

def test_long_read_dna_variant_calling_deepvariant_local():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_dna_bam_file = get_data_path(name='bam/sample001tumor_minimap2_mdtagged_sorted.bam')
    long_read_tumor_dna_bam_bai_file = get_data_path(name='bam/sample001tumor_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_dna_variant_calling_deepvariant'
    work_dir = temp_dir + '/work/test_long_read_dna_variant_calling_deepvariant'
    output_dir = temp_dir + '/outputs/test_long_read_dna_variant_calling_deepvariant'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample001tumor'],
        'bam_file': [long_read_tumor_dna_bam_file],
        'bam_bai_file': [long_read_tumor_dna_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--containerization', 'docker',
        '--deepvariant_input_path', '/Users/leework/Documents/Research/projects/project_nexus/',
        '--deepvariant_output_path', '/var/folders/',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='long_read_dna_variant_calling_deepvariant.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


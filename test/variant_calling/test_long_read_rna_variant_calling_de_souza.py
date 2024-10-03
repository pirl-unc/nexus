import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_long_read_rna_variant_calling_de_souza_github():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_fastq_file = get_data_path(name='fastq/sample200tumor_long_read_rna.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa')
    reference_genome_fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa.fai')
    reference_genome_fasta_dict_file = get_data_path(name='fasta/hg38_chr17_1-8M.dict')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_rna_variant_calling_de_souza'
    work_dir = temp_dir + '/work/test_long_read_rna_variant_calling_de_souza'
    output_dir = temp_dir + '/outputs/test_long_read_rna_variant_calling_de_souza'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample200tumor'],
        'fastq_file': [tumor_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--reference_genome_fasta_dict_file', reference_genome_fasta_dict_file,
        '--deepvariant_containerization', 'docker',
        '--deepvariant_input_path', '/home/runner/work/nexus/nexus/',
        '--deepvariant_output_path', '/tmp/',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='long_read_rna_variant_calling_de_souza.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_long_read_rna_variant_calling_de_souza_local():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_fastq_file = get_data_path(name='fastq/sample200tumor_long_read_rna.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa')
    reference_genome_fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa.fai')
    reference_genome_fasta_dict_file = get_data_path(name='fasta/hg38_chr17_1-8M.dict')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_rna_variant_calling_de_souza'
    work_dir = temp_dir + '/work/test_long_read_rna_variant_calling_de_souza'
    output_dir = temp_dir + '/outputs/test_long_read_rna_variant_calling_de_souza'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample200tumor'],
        'fastq_file': [tumor_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--reference_genome_fasta_dict_file', reference_genome_fasta_dict_file,
        '--deepvariant_containerization', 'docker',
        '--deepvariant_input_path', '/Users/leework/Documents/Research/projects/project_nexus/',
        '--deepvariant_output_path', '/var/folders/',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='long_read_rna_variant_calling_de_souza.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

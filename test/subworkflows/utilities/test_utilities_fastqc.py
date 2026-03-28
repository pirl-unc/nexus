import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_utilities_fastqc_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fastq_file = get_data_path(name='fastq/nexus-dna-001-tumor_long_read.fastq.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_utilities_fastqc_1'
    work_dir = temp_dir + '/work/test_utilities_fastqc_1'
    output_dir = temp_dir + '/outputs/test_utilities_fastqc_1'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor_long_reads'],
        'fastq_file_1': [fastq_file]
    }).to_csv(intermediate_dir + "/samples1.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples1.tsv',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='utilities_fastqc.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_utilities_fastqc_2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    dna_fastq_file_1 = get_data_path(name='fastq/nexus-dna-001-tumor_paired-end_read_r1.fastq.gz')
    dna_fastq_file_2 = get_data_path(name='fastq/nexus-dna-001-tumor_paired-end_read_r2.fastq.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_utilities_fastqc_2'
    work_dir = temp_dir + '/work/test_utilities_fastqc_2'
    output_dir = temp_dir + '/outputs/test_utilities_fastqc_2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor_short_reads'],
        'fastq_file_1': [dna_fastq_file_1],
        'fastq_file_2': [dna_fastq_file_2]
    }).to_csv(intermediate_dir + "/samples2.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples2.tsv',
        '--read_type', 'paired-end',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='utilities_fastqc.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

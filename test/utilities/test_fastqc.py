import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_fastqc_single_end_reads():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_dna_fastq_file = get_data_path(name='fastq/hg38_tp53_tumor_long_read_dna.fastq.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_fastqc'
    work_dir = temp_dir + '/work/test_fastqc'
    output_dir = temp_dir + '/outputs/test_fastqc'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor-single-end'],
        'fastq_file_1': [long_read_tumor_dna_fastq_file]
    }).to_csv(intermediate_dir + "/samples1.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples1.tsv',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='fastqc.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_fastqc_paired_end_reads():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    dna_fastq_file_1 = get_data_path(name='fastq/hg38_tp53_tumor_paired-end_read_dna.r1.fastq.gz')
    dna_fastq_file_2 = get_data_path(name='fastq/hg38_tp53_tumor_paired-end_read_dna.r2.fastq.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_fastqc'
    work_dir = temp_dir + '/work/test_fastqc'
    output_dir = temp_dir + '/outputs/test_fastqc'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor-paired-end'],
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
    run_workflow(workflow='fastqc.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

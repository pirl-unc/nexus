import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path, get_alias_path


def test_docker_paired_end_read_rna_hla_typing_seq2hla():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fastq_file_1 = get_data_path(name='fastq/sample300normal_paired-end_read_rna.r1.fastq.gz')
    fastq_file_2 = get_data_path(name='fastq/sample300normal_paired-end_read_rna.r2.fastq.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_paired_end_rna_hla_profiling_seq2hla'
    work_dir = temp_dir + '/work/test_paired_end_rna_hla_profiling_seq2hla'
    output_dir = temp_dir + '/outputs/test_paired_end_rna_hla_profiling_seq2hla'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample300normal'],
        'fastq_file_1': [fastq_file_1],
        'fastq_file_2': [fastq_file_2]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='paired-end_read_rna_hla_typing_seq2hla.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


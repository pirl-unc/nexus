import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_assembly_spades():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fastq_file_1 = get_data_path(name='fastq/nexus-rna-001-tumor_paired-end_read_r1.fastq.gz')
    fastq_file_2 = get_data_path(name='fastq/nexus-rna-001-tumor_paired-end_read_r2.fastq.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_assembly_spades'
    work_dir = temp_dir + '/work/test_assembly_spades'
    output_dir = temp_dir + '/outputs/test_assembly_spades'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-001-tumor'],
        'fastq_file_1': [fastq_file_1],
        'fastq_file_2': [fastq_file_2]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='assembly_spades.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_hla_typing_hlaminer_sr_dna():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fastq_file_1 = get_data_path(name='fastq/nexus-dna-104-normal_paired-end_read_r1.fastq.gz')
    fastq_file_2 = get_data_path(name='fastq/nexus-dna-104-normal_paired-end_read_r2.fastq.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_hla_typing_hlaminer-sr-dna'
    work_dir = temp_dir + '/work/test_hla_typing_hlaminer-sr-dna'
    output_dir = temp_dir + '/outputs/test_hla_typing_hlaminer-sr-dna'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-104-normal'],
        'fastq_file_1': [fastq_file_1],
        'fastq_file_2': [fastq_file_2]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--params_bwa_aln', '"-e 0 -o 0"',
        '--params_bwa_sampe', '"-o 1000"',
        '--params_hlaminer', '"-s 500"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='hla_typing_hlaminer-sr-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

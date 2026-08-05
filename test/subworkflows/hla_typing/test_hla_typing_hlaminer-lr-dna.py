import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_hla_typing_hlaminer_lr_dna():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fastq_file = get_data_path(name='fastq/nexus-dna-104-normal_long_read.fastq.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_hla_typing_hlaminer-lr-dna'
    work_dir = temp_dir + '/work/test_hla_typing_hlaminer-lr-dna'
    output_dir = temp_dir + '/outputs/test_hla_typing_hlaminer-lr-dna'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-104-normal'],
        'fastq_file': [fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--params_minimap2', '"-ax map-ont --secondary=no"',
        '--params_hlaminer', '"-s 500 -q 1 -i 1"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='hla_typing_hlaminer-lr-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

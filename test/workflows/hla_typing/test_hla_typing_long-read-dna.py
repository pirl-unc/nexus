import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


# FuFiHLA, SpecHLA (long-read), SpecImmune (nexus-dna-104-normal)
def test_hla_typing_longread_dna():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fastq_file = get_data_path(name='fastq/nexus-dna-104-normal_long_read.fastq.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/hla_typing_long-read-dna/params_1.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_hla_typing_longread_dna'
    work_dir = temp_dir + '/work/test_hla_typing_longread_dna'
    output_dir = temp_dir + '/outputs/test_hla_typing_longread_dna'
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

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir

    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='hla_typing_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

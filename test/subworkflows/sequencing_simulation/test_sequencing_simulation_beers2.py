import os
import pandas as pd
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_sequencing_simulation_beers2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    config_yaml_file = get_data_path(name='inputs/subworkflows/beers2/input_data/test_sample/beers2.config.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_sequencing_simulation_beers2'
    work_dir = temp_dir + '/work/test_sequencing_simulation_beers2'
    output_dir = temp_dir + '/outputs/test_sequencing_simulation_beers2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['test_sample'],
        'config_file': [config_yaml_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='sequencing_simulation_beers2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

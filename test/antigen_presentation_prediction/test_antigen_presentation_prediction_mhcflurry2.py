import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_antigen_presentation_prediction_mhcflurry2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    csv_file = get_data_path(name='csv/sample900_mhcflurry2_input.csv')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_antigen_presentation_prediction_mhcflurry2'
    work_dir = temp_dir + '/work/test_antigen_presentation_prediction_mhcflurry2'
    output_dir = temp_dir + '/outputs/test_antigen_presentation_prediction_mhcflurry2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample900'],
        'mhcflurry2_input_csv_file': [csv_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='antigen_presentation_prediction_mhcflurry2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


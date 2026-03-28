import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_utilities_pbccs():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='inputs/subworkflows/pbccs/sample_pbsim3.bam')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_utilities_pbccs_1'
    work_dir = temp_dir + '/work/test_utilities_pbccs_1'
    output_dir = temp_dir + '/outputs/test_utilities_pbccs_1'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['test'],
        'subreads_bam_file': [bam_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='utilities_pbccs.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

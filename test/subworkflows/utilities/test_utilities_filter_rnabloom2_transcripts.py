import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_utilities_filter_rnabloom2_transcripts():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    rnabloom2_output_dir = get_data_path(name='inputs/subworkflows/utilities_filter-rnabloom2-transcripts/nexus-rna-002-tumor_rnabloom2_output/')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_utilities_filter_rnabloom2_transcripts'
    work_dir = temp_dir + '/work/test_utilities_filter_rnabloom2_transcripts'
    output_dir = temp_dir + '/outputs/test_utilities_filter_rnabloom2_transcripts'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-002-tumor'],
        'rnabloom2_output_dir': [rnabloom2_output_dir]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='utilities_filter-rnabloom2-transcripts.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

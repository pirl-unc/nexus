import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_assembly_isonform_clustered():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fastq_file = get_data_path(name='inputs/subworkflows/isonform_clustered/reads.fastq.gz')
    tsv_file = get_data_path(name='inputs/subworkflows/isonform_clustered/groups.tsv')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_assembly_isonform'
    work_dir = temp_dir + '/work/test_assembly_isonform'
    output_dir = temp_dir + '/outputs/test_assembly_isonform'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample'],
        'fastq_file': [fastq_file],
        'tsv_file': [tsv_file],
        'cluster_method': ['isonclust3']
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='assembly_isonform_clustered.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

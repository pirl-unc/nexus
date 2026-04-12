import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_antigen_prediction_netmhcpan4():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fasta_file = get_data_path(name='fasta/sample_netmhcpan_input.fasta')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_antigen_prediction_netmhcpan4'
    work_dir = temp_dir + '/work/test_antigen_prediction_netmhcpan4'
    output_dir = temp_dir + '/outputs/test_antigen_prediction_netmhcpan4'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample'],
        'fasta_file': [fasta_file],
        'mhc_alleles': ['H-2-Kb,H-2-Db']
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--netmhcpan4_home_dir', '/Users/leework/Documents/Software/netmhc-bundle-master/netMHCpan-4.1/Darwin_x86_64/',
        '--netmhcpan4', '/Users/leework/Documents/Software/netmhc-bundle-master/netMHCpan-4.1/Darwin_x86_64/bin/netMHCpan',
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='antigen_prediction_netmhcpan4.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

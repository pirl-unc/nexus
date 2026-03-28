import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_sequencing_simulation_pbsim3_dna_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fasta_file = get_data_path(name='fasta/sample_genome.fa')
    pbsim3_model_file = get_data_path(name='indices/pbsim3/ERRHMM-SEQUEL.model')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_sequencing_simulation_pbsim3_dna_1'
    work_dir = temp_dir + '/work/test_sequencing_simulation_pbsim3_dna_1'
    output_dir = temp_dir + '/outputs/test_sequencing_simulation_pbsim3_dna_1'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['test_dna_1'],
        'fasta_file': [fasta_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--pbsim3_model_file', pbsim3_model_file,
        '--params_pbsim3_mode', '"errhmm"',
        '--params_pbsim3', '"--depth 1 --length-mean 10000 --length-sd 3000 --pass-num 5 --strategy wgs"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='sequencing_simulation_pbsim3-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_sequencing_simulation_pbsim3_dna_2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fasta_file = get_data_path(name='fasta/sample_genome.fa')
    pbsim3_model_file = get_data_path(name='indices/pbsim3/QSHMM-ONT-HQ.model')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_sequencing_simulation_pbsim3_dna_2'
    work_dir = temp_dir + '/work/test_sequencing_simulation_pbsim3_dna_2'
    output_dir = temp_dir + '/outputs/test_sequencing_simulation_pbsim3_dna_2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['test_dna_2'],
        'fasta_file': [fasta_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--pbsim3_model_file', pbsim3_model_file,
        '--params_pbsim3_mode', '"qshmm"',
        '--params_pbsim3', '"--depth 1 --length-mean 10000 --length-sd 3000 --strategy wgs"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='sequencing_simulation_pbsim3-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


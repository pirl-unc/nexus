import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_sequencing_simulation_nanosim_genome():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fasta_file = get_data_path(name='fasta/sample_genome.fa')
    model_prefix = get_data_path(name='indices/nanosim/human_giab_hg002_sub1M_kitv14_dorado_v3.2.1/training')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_sequencing_simulation_nanosim_genome_1'
    work_dir = temp_dir + '/work/test_sequencing_simulation_nanosim_genome_1'
    output_dir = temp_dir + '/outputs/test_sequencing_simulation_nanosim_genome_1'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['test1_b1'],
        'fasta_file': [fasta_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--model_prefix', model_prefix,
        '--params_nanosim', '"--coverage 1 --homopolymer --KmerBias 5"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='sequencing_simulation_nanosim-genome.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

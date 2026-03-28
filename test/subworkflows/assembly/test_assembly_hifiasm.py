import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_assembly_hifiasm():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fastx_file = get_data_path(name='fasta/hifiasm_demo_chr11-2M.fa.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_assembly_hifiasm'
    work_dir = temp_dir + '/work/test_assembly_hifiasm'
    output_dir = temp_dir + '/outputs/test_assembly_hifiasm'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample'],
        'fastq_file': [fastx_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--params_hifiasm', '"--primary"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='assembly_hifiasm.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

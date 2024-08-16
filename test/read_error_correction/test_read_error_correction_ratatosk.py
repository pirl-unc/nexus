import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_read_error_correction_ratatosk():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_dna_fastq_file = get_data_path(name='fastq/sample001tumor_long_read_dna.fastq.gz')
    short_read_tumor_dna_fastq_r1_file = get_data_path(name='fastq/sample100tumor_paired-end_read_dna.r1.fastq.gz')
    short_read_tumor_dna_fastq_r2_file = get_data_path(name='fastq/sample100tumor_paired-end_read_dna.r2.fastq.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_read_error_correction_ratatosk'
    work_dir = temp_dir + '/work/test_read_error_correction_ratatosk'
    output_dir = temp_dir + '/outputs/test_read_error_correction_ratatosk'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample001tumor'],
        'long_read_fastq_file': [long_read_tumor_dna_fastq_file],
        'short_read_r1_fastq_file': [short_read_tumor_dna_fastq_r1_file],
        'short_read_r2_fastq_file': [short_read_tumor_dna_fastq_r2_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='long_read_error_correction_ratatosk.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


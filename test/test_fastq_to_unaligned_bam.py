import pandas as pd
import os
from .data import get_data_path, get_alias_path
from nexuslib.main import run_workflow


def test_fastq_to_unaligned_bam():
    picard_jar_path = get_alias_path(executable_name='picard.jar')
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    long_read_tumor_dna_fastq_file = get_data_path(name='hg38_tp53_tumor_long_read_dna.fastq.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_fastq_to_unaligned_bam'
    work_dir = temp_dir + '/work/test_fastq_to_unaligned_bam'
    output_dir = temp_dir + '/outputs/test_fastq_to_unaligned_bam'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor'],
        'fastq_file_1': [long_read_tumor_dna_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--read_type', 'single-end',
        '--picard', picard_jar_path,
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='fastq_to_unaligned_bam.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_sequencing_coverage():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_dna_bam_file = get_data_path(name='bam/hg38_tp53_tumor_long_read_dna.bam')
    long_read_tumor_dna_bam_bai_file = get_data_path(name='bam/hg38_tp53_tumor_long_read_dna.bam.bai')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_sequencing_coverage'
    work_dir = temp_dir + '/work/test_sequencing_coverage'
    output_dir = temp_dir + '/outputs/test_sequencing_coverage'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor'],
        'bam_file': [long_read_tumor_dna_bam_file],
        'bam_bai_file': [long_read_tumor_dna_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='sequencing_coverage.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

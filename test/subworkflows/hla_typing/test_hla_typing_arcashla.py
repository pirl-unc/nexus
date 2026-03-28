import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_hla_typing_arcashla():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-rna-104-normal-paired-end-read_star_outputs/nexus-rna-104-normal-paired-end-read_star_Aligned.sortedByCoord.out.bam')
    bam_bai_file = get_data_path(name='bam/nexus-rna-104-normal-paired-end-read_star_outputs/nexus-rna-104-normal-paired-end-read_star_Aligned.sortedByCoord.out.bam.bai')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_hla_typing_arcashla'
    work_dir = temp_dir + '/work/test_hla_typing_arcashla'
    output_dir = temp_dir + '/outputs/test_hla_typing_arcashla'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-104-normal'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='hla_typing_arcashla.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

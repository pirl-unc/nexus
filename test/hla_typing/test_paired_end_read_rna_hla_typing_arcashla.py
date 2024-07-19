import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_paired_end_read_rna_hla_typing_arcashla():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/sample001_star_Aligned.sortedByCoord.out.bam')
    bam_bai_file = get_data_path(name='bam/sample001_star_Aligned.sortedByCoord.out.bam.bai')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_paired_end_rna_hla_profiling_arcashla'
    work_dir = temp_dir + '/work/test_paired_end_rna_hla_profiling_arcashla'
    output_dir = temp_dir + '/outputs/test_paired_end_rna_hla_profiling_arcashla'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample001'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='paired-end_read_rna_hla_typing_arcashla.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

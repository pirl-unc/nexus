import pandas as pd
import os
from .data import get_data_path
from nexuslib.main import run_workflow


def test_novel_isoform_discovery_isoseq():
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    long_read_tumor_rna_unaligned_bam_file = get_data_path(name='hg38_tp53_tumor_long_read_rna_unaligned.bam')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_novel_isoform_discovery_isoseq'
    work_dir = temp_dir + '/work/test_novel_isoform_discovery_isoseq'
    output_dir = temp_dir + '/outputs/test_novel_isoform_discovery_isoseq'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor'],
        'bam_file': [long_read_tumor_rna_unaligned_bam_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--isoseq', 'isoseq',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='novel_isoform_discovery_isoseq.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


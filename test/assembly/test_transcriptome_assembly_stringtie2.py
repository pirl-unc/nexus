import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_transcriptome_assembly_stringtie2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/sample200tumor_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/sample200tumor_minimap2_mdtagged_sorted.bam.bai')
    gtf_file = get_data_path(name='gtf/gencode_v41_tp53_annotation.gtf')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_transcriptome_assembly_stringtie2'
    work_dir = temp_dir + '/work/test_transcriptome_assembly_stringtie2'
    output_dir = temp_dir + '/outputs/test_transcriptome_assembly_stringtie2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample200tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--gtf_file', gtf_file,
        '--output_dir', output_dir
    ]
    run_workflow(workflow='transcriptome_assembly_stringtie2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

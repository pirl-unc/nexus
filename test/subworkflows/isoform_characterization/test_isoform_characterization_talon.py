import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_isoform_characterization_talon():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_bam_file = get_data_path(name='bam/nexus-rna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    long_read_tumor_bam_bai_file = get_data_path(name='bam/nexus-rna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_isoform_characterization_talon'
    work_dir = temp_dir + '/work/test_isoform_characterization_talon'
    output_dir = temp_dir + '/outputs/test_isoform_characterization_talon'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-002-tumor'],
        'bam_file': [long_read_tumor_bam_file],
        'bam_bai_file': [long_read_tumor_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--output_dir', output_dir
    ]
    run_workflow(workflow='isoform_characterization_talon.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

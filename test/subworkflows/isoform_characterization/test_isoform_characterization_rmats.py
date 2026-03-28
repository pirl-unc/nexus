import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_isoform_characterization_rmats():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-rna-103-tumor-paired-end-read_star_outputs/nexus-rna-103-tumor-paired-end-read_star_Aligned.sortedByCoord.out.bam')
    bam_bai_file = get_data_path(name='bam/nexus-rna-103-tumor-paired-end-read_star_outputs/nexus-rna-103-tumor-paired-end-read_star_Aligned.sortedByCoord.out.bam.bai')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_isoform_characterization_rmats'
    work_dir = temp_dir + '/work/test_isoform_characterization_rmats'
    output_dir = temp_dir + '/outputs/test_isoform_characterization_rmats'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-103-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--params_rmats', '"-t paired --libType fr-unstranded --readLength 151 --mil 10 --novelSS --statoff"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='isoform_characterization_rmats.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

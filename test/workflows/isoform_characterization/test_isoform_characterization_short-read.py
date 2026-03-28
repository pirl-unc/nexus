import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


# rMATS (nexus-rna-103-tumor)
def test_isoform_characterization_short_read():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-rna-103-tumor-paired-end-read_star_outputs/nexus-rna-103-tumor-paired-end-read_star_Aligned.sortedByCoord.out.bam')
    bam_bai_file = get_data_path(name='bam/nexus-rna-103-tumor-paired-end-read_star_outputs/nexus-rna-103-tumor-paired-end-read_star_Aligned.sortedByCoord.out.bam.bai')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/isoform_characterization_short-read/params_1.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_isoform_characterization_shortread_1'
    work_dir = temp_dir + '/work/test_isoform_characterization_shortread_1'
    output_dir = temp_dir + '/outputs/test_isoform_characterization_shortread_1'
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

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genes_gtf_file'] = reference_genes_gtf_file

    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='isoform_characterization_short-read.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_assembly_long_read_rna():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fastq_file = get_data_path(name='inputs/workflows/assembly_long-read-rna/reads.fastq.gz')
    tsv_file = get_data_path(name='inputs/workflows/assembly_long-read-rna/groups.tsv')
    bam_file = get_data_path(name='bam/nexus-rna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-rna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/assembly_long-read-rna/params.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_assembly_long_read_rna'
    work_dir = temp_dir + '/work/test_assembly_long_read_rna'
    output_dir = temp_dir + '/outputs/test_assembly_long_read_rna'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample'],
        'fastq_file': [fastq_file],
        'tsv_file': [tsv_file],
        'cluster_method': ['isonclust3'],
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
    run_workflow(workflow='assembly_long-read-rna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

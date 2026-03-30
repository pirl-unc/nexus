import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


# All tools + SQANTI3-GTF chaining (nexus-rna-003-tumor)
def test_isoform_characterization_long_read():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-rna-103-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-rna-103-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    fastq_file = get_data_path(name='fastq/nexus-rna-103-tumor_long_read.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/isoform_characterization_long-read/params.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_isoform_characterization_longread'
    work_dir = temp_dir + '/work/test_isoform_characterization_longread'
    output_dir = temp_dir + '/outputs/test_isoform_characterization_longread'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-103-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file],
        'fastq_file': [fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['reference_genes_gtf_file'] = reference_genes_gtf_file

    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='isoform_characterization_long-read.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

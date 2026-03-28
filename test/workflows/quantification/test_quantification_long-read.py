import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


# Bambu, LIQA (BAM-based — nexus-rna-001-tumor)
def test_quantification_long_read_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-rna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-rna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    fastq_file = get_data_path(name='fastq/nexus-rna-002-tumor_long_read.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    reference_transcripts_fasta_file = get_data_path(name='fasta/gencode.v45.transcripts.tp53.fa')
    params_yaml_file = get_data_path(name='inputs/workflows/quantification_long-read/params_1.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_quantification_longread_1'
    work_dir = temp_dir + '/work/test_quantification_longread_1'
    output_dir = temp_dir + '/outputs/test_quantification_longread_1'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-002-tumor'],
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
    params['reference_transcripts_fasta_file'] = reference_transcripts_fasta_file

    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='quantification_long-read.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# Kallisto, Oarfish, TranSigner (FASTQ-based — nexus-rna-002-tumor)
def test_quantification_long_read_2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-rna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-rna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    fastq_file = get_data_path(name='fastq/nexus-rna-002-tumor_long_read.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    reference_transcripts_fasta_file = get_data_path(name='fasta/gencode.v45.transcripts.tp53.fa')
    params_yaml_file = get_data_path(name='inputs/workflows/quantification_long-read/params_2.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_quantification_longread_2'
    work_dir = temp_dir + '/work/test_quantification_longread_2'
    output_dir = temp_dir + '/outputs/test_quantification_longread_2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-002-tumor'],
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
    params['reference_transcripts_fasta_file'] = reference_transcripts_fasta_file

    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='quantification_long-read.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

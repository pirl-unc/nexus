import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_novel_isoform_discovery_sqanti3_gtf_mode_flair():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_rna_gtf_file = get_data_path(name='gtf/sample200tumor_flair.isoforms.gtf')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa')
    reference_genome_fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa.fai')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode_v41_tp53_annotation.gtf')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_novel_isoform_discovery_sqanti3_gtf_mode'
    work_dir = temp_dir + '/work/test_novel_isoform_discovery_sqanti3_gtf_mode'
    output_dir = temp_dir + '/outputs/test_novel_isoform_discovery_sqanti3_gtf_mode'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample200tumor_flair'],
        'gtf_file': [long_read_tumor_rna_gtf_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--params_sqanti3_qc', '"--report skip"',
        '--params_sqanti3_filter', '"ml"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='novel_isoform_discovery_sqanti3_gtf_mode.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_novel_isoform_discovery_sqanti3_gtf_mode_isoquant():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_rna_gtf_file = get_data_path(name='gtf/sample200tumor_isoquant.transcript_models.gtf')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa')
    reference_genome_fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa.fai')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode_v41_tp53_annotation.gtf')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_novel_isoform_discovery_sqanti3_gtf_mode'
    work_dir = temp_dir + '/work/test_novel_isoform_discovery_sqanti3_gtf_mode'
    output_dir = temp_dir + '/outputs/test_novel_isoform_discovery_sqanti3_gtf_mode'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample200tumor_isoquant'],
        'gtf_file': [long_read_tumor_rna_gtf_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--params_sqanti3_qc', '"--report skip"',
        '--params_sqanti3_filter', '"ml"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='novel_isoform_discovery_sqanti3_gtf_mode.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

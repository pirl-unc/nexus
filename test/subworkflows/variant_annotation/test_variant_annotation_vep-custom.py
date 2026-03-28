import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_annotation_vep_custom_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    vcf_file = get_data_path(name='vcf/nexus-dna-001-tumor_gatk4-mutect2.vcf')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_annotation_vep_custom_1'
    work_dir = temp_dir + '/work/test_variant_annotation_vep_custom_1'
    output_dir = temp_dir + '/outputs/test_variant_annotation_vep_custom_1'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'vcf_file': [vcf_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--vep_dir', '/Users/leework/.vep',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--reference_genes_gtf_source', 'GENCODE',
        '--params_vep', '"--species homo_sapiens --offline --cache --assembly GRCh38"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_annotation_vep-custom.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_variant_annotation_vep_custom_2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    vcf_file = get_data_path(name='vcf/nexus-dna-001-tumor_strelka2_snvs.vcf')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_annotation_vep_custom_2'
    work_dir = temp_dir + '/work/test_variant_annotation_vep_custom_2'
    output_dir = temp_dir + '/outputs/test_variant_annotation_vep_custom_2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'vcf_file': [vcf_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--vep_dir', '/Users/leework/.vep',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--reference_genes_gtf_source', 'GENCODE',
        '--params_vep', '"--species homo_sapiens --offline --cache --assembly GRCh38"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_annotation_vep-custom.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_variant_annotation_vep_custom_3():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    vcf_file = get_data_path(name='vcf/nexus-dna-002-tumor_strelka2_indels.vcf')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_annotation_vep_custom_3'
    work_dir = temp_dir + '/work/test_variant_annotation_vep_custom_3'
    output_dir = temp_dir + '/outputs/test_variant_annotation_vep_custom_3'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-002-tumor'],
        'vcf_file': [vcf_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--vep_dir', '/Users/leework/.vep',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--reference_genes_gtf_source', 'GENCODE',
        '--params_vep', '"--species homo_sapiens --offline --cache --assembly GRCh38"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_annotation_vep-custom.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


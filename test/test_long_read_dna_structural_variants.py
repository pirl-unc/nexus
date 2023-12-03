import pandas as pd
import os
from .data import get_data_path
from nexuslib.main import run_workflow


def test_long_read_dna_structural_variants_local():
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    long_read_tumor_dna_bam_file = get_data_path(name='hg38_tp53_tumor_long_read_dna_minimap2_mdtagged_sorted.bam')
    long_read_tumor_dna_bam_bai_file = get_data_path(name='hg38_tp53_tumor_long_read_dna_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='hg38_chr17_1-8000000.fa')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_dna_structural_variants'
    work_dir = temp_dir + '/work/test_long_read_dna_structural_variants'
    output_dir = temp_dir + '/outputs/test_long_read_dna_structural_variants'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor'],
        'bam_file': [long_read_tumor_dna_bam_file],
        'bam_bai_file': [long_read_tumor_dna_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--tools_list', 'sniffles2,svim,cutesv',    # pbsv is only available on linux
        '--output_dir', output_dir
    ]
    run_workflow(workflow='long_read_dna_structural_variants.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_long_read_dna_structural_variants_github():
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    long_read_tumor_dna_bam_file = get_data_path(name='hg38_tp53_tumor_long_read_dna_minimap2_mdtagged_sorted.bam')
    long_read_tumor_dna_bam_bai_file = get_data_path(name='hg38_tp53_tumor_long_read_dna_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='hg38_chr17_1-8000000.fa')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_dna_structural_variants'
    work_dir = temp_dir + '/work/test_long_read_dna_structural_variants'
    output_dir = temp_dir + '/outputs/test_long_read_dna_structural_variants'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor'],
        'bam_file': [long_read_tumor_dna_bam_file],
        'bam_bai_file': [long_read_tumor_dna_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--tools_list', 'sniffles2,svim,pbsv,cutesv',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='long_read_dna_structural_variants.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

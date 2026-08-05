import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_haplotagging_longcallr_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_rna_bam_file = get_data_path(name='bam/nexus-rna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    long_read_tumor_rna_bam_bai_file = get_data_path(name='bam/nexus-rna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_longcallr_1'
    work_dir = temp_dir + '/work/test_haplotagging_longcallr_1'
    output_dir = temp_dir + '/outputs/test_haplotagging_longcallr_1'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-001-tumor'],
        'bam_file': [long_read_tumor_rna_bam_file],
        'bam_bai_file': [long_read_tumor_rna_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--preset', 'hifi-masseq',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='haplotagging_longcallr.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_haplotagging_longcallr_2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='inputs/subworkflows/longcallr/demo.bam')
    bai_file = get_data_path(name='inputs/subworkflows/longcallr/demo.bam.bai')
    fasta_file = get_data_path(name='inputs/subworkflows/longcallr/chr20.fa')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr20.gtf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_longcallr_2'
    work_dir = temp_dir + '/work/test_haplotagging_longcallr_2'
    output_dir = temp_dir + '/outputs/test_haplotagging_longcallr_2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['demo'],
        'bam_file': [bam_file],
        'bam_bai_file': [bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', fasta_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--preset', 'hifi-masseq',
        '--haplotag_output', '"both"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='haplotagging_longcallr.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

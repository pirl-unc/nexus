import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_calling_isolaser():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='inputs/subworkflows/isolaser/dlpfc_1.bam')
    bam_bai_file = get_data_path(name='inputs/subworkflows/isolaser/dlpfc_1.bam.bai')
    fastq_file = get_data_path(name='inputs/subworkflows/isolaser/dlpfc_1.fq.gz')
    gtf_file = get_data_path(name='inputs/subworkflows/isolaser/gencode.v32.chr21.gtf.gz')
    reference_genome_fasta_file = get_data_path(name='inputs/subworkflows/isolaser/GRCh38.chr21.fa')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_isolaser'
    work_dir = temp_dir + '/work/test_variant_calling_isolaser'
    output_dir = temp_dir + '/outputs/test_variant_calling_isolaser'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['dlpfc_1'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file],
        'fastq_file': [fastq_file],
        'gtf_file': [gtf_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_isolaser.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_alignment_ultra():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_rna_fastq_file_1 = get_data_path(name='fastq/nexus-rna-001-tumor_long_read.fastq.gz')
    long_read_rna_fastq_file_2 = get_data_path(name='fastq/nexus-rna-002-tumor_long_read.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_alignment_ultra'
    work_dir = temp_dir + '/work/test_alignment_ultra/'
    output_dir = temp_dir + '/outputs/test_alignment_ultra/'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-001-tumor', 'nexus-rna-002-tumor'],
        'fastq_file': [long_read_rna_fastq_file_1, long_read_rna_fastq_file_2]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    # Run uLTRA
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--output_dir', output_dir
    ]
    run_workflow(workflow='alignment_ultra.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

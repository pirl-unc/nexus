import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_long_read_rna_alignment_ultra():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_rna_fastq_file = get_data_path(name='fastq/sample200tumor_long_read_rna.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa')
    reference_genome_fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa.fai')
    ultra_index = get_data_path(name='indices/ultra/hg38_chr17/')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_rna_alignment_ultra'
    work_dir = temp_dir + '/work/test_long_read_rna_alignment_ultra/'
    output_dir = temp_dir + '/outputs/test_long_read_rna_alignment_ultra/'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor'],
        'fastq_file': [long_read_tumor_rna_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    # Run uLTRA
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--ultra_index', ultra_index,
        '--output_dir', output_dir
    ]
    run_workflow(workflow='long_read_rna_alignment_ultra.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

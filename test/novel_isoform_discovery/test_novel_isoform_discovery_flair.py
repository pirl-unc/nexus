import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_novel_isoform_discovery_flair():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_rna_fastq_file = get_data_path(name='fastq/hg38_tp53_tumor_long_read_rna.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa')
    gtf_file = get_data_path(name='gtf/gencode_v41_tp53_annotation.gtf')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_novel_isoform_discovery_flair'
    work_dir = temp_dir + '/work/test_novel_isoform_discovery_flair'
    output_dir = temp_dir + '/outputs/test_novel_isoform_discovery_flair'
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
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--gtf_file', gtf_file,
        '--output_dir', output_dir
    ]
    run_workflow(workflow='novel_isoform_discovery_flair.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

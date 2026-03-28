import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_quantification_kallisto_lr():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fastq_file = get_data_path(name='fastq/nexus-rna-002-tumor_long_read.fastq.gz')
    reference_transcripts_fasta_file = get_data_path(name='fasta/gencode.v45.transcripts.tp53.fa')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_quantification_kallisto_lr'
    work_dir = temp_dir + '/work/test_quantification_kallisto_lr'
    output_dir = temp_dir + '/outputs/test_quantification_kallisto_lr'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-002-tumor'],
        'fastq_file': [fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_transcripts_fasta_file', reference_transcripts_fasta_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--output_dir', output_dir
    ]
    run_workflow(workflow='quantification_kallisto-lr.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

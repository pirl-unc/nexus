import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_quantification_salmon_fastq():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    paired_end_read_tumor_fastq_file_1 = get_data_path(name='fastq/nexus-rna-002-tumor_paired-end_read_r1.fastq.gz')
    paired_end_read_tumor_fastq_file_2 = get_data_path(name='fastq/nexus-rna-002-tumor_paired-end_read_r2.fastq.gz')
    reference_transcripts_fasta_file = get_data_path(name='fasta/gencode.v45.transcripts.tp53.fa')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_quantification_salmon_fastq'
    work_dir = temp_dir + '/work/test_quantification_salmon_fastq'
    output_dir = temp_dir + '/outputs/test_quantification_salmon_fastq'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-002-tumor'],
        'fastq_file_1': [paired_end_read_tumor_fastq_file_1],
        'fastq_file_2': [paired_end_read_tumor_fastq_file_2]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_transcripts_fasta_file', reference_transcripts_fasta_file,
        '--params_salmon_index', '"--gencode "',
        '--params_salmon_quant', '"--libType IU --seqBias --gcBias --posBias"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='quantification_salmon-fastq.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

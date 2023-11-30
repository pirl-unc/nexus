import pandas as pd
import os
from .data import get_data_path
from nexuslib.main import run_workflow


def test_long_read_alignment_minimap2_dna():
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    long_read_tumor_dna_fastq_file = get_data_path(name='hg38_tp53_tumor_long_read_dna.fastq.gz')
    long_read_normal_dna_fastq_file = get_data_path(name='hg38_tp53_normal_long_read_dna.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='hg38_chr17_1-8000000.fa')
    temp_dir = os.getcwd() + '/tmp/'
    intermediate_dir = temp_dir + 'intermediate/test_long_read_alignment_minimap2_dna/'
    work_dir = temp_dir + 'work/test_long_read_alignment_minimap2_dna/'
    output_dir = temp_dir + 'outputs/test_long_read_alignment_minimap2_dna/'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor', 'normal'],
        'fastq_file': [long_read_tumor_dna_fastq_file, long_read_normal_dna_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--minimap2', 'minimap2',
        '--minimap2_params', '"-ax map-hifi --cs --eqx -Y -L "',
        '--samtools', 'samtools',
        '--platform_tag', 'pacbio',
        '--output_dir', output_dir,
        '-w', work_dir
    ]
    run_workflow(workflow='long_read_alignment_minimap2.nf',
                 workflow_args=workflow_args)


def test_long_read_alignment_minimap2_rna():
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    long_read_tumor_rna_fastq_file = get_data_path(name='hg38_tp53_tumor_long_read_rna.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='hg38_chr17_1-8000000.fa')
    temp_dir = os.getcwd() + '/tmp/'
    intermediate_dir = temp_dir + 'intermediate/test_alignment_minimap2_rna/'
    work_dir = temp_dir + 'work/test_alignment_minimap2_rna/'
    output_dir = temp_dir + 'outputs/test_alignment_minimap2_rna/'
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
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--minimap2', 'minimap2',
        '--minimap2_params', "'-ax splice:hq -uf --cs --eqx -Y -L '",
        '--samtools', 'samtools',
        '--platform_tag', 'pacbio',
        '--output_dir', output_dir,
        '-w', work_dir
    ]
    run_workflow(workflow='long_read_alignment_minimap2.nf',
                 workflow_args=workflow_args)

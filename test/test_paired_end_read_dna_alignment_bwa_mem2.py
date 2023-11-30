import pandas as pd
import os
from .data import get_data_path
from nexuslib.main import run_workflow


def test_paired_end_read_dna_alignment_bwa_mem2():
    abra2_jar_path = os.environ.get('abra2')
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    paired_end_read_tumor_dna_r1_fastq_file = get_data_path(name='hg38_tp53_tumor_paired_end_dna.r1.fastq.gz')
    paired_end_read_tumor_dna_r2_fastq_file = get_data_path(name='hg38_tp53_tumor_paired_end_dna.r2.fastq.gz')
    paired_end_read_normal_dna_r1_fastq_file = get_data_path(name='hg38_tp53_normal_paired_end_dna.r1.fastq.gz')
    paired_end_read_normal_dna_r2_fastq_file = get_data_path(name='hg38_tp53_normal_paired_end_dna.r2.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='hg38_chr17_1-8000000.fa')
    abra2_targets_bed_file = get_data_path(name='gencode-v41-annotation-abra2-exon-targets.bed')
    known_sites_vcf_file = get_data_path(name='known_sites_hg38.vcf')
    temp_dir = os.getcwd() + '/tmp/'
    intermediate_dir = temp_dir + 'intermediate/test_paired_end_read_alignment_bwa_mem2/'
    work_dir = temp_dir + 'work/test_paired_end_read_alignment_bwa_mem2/'
    output_dir = temp_dir + 'outputs/test_paired_end_read_alignment_bwa_mem2/'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor', 'normal'],
        'fastq_file_1': [paired_end_read_tumor_dna_r1_fastq_file, paired_end_read_normal_dna_r1_fastq_file],
        'fastq_file_2': [paired_end_read_tumor_dna_r2_fastq_file, paired_end_read_normal_dna_r2_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--bwa_mem2', 'bwa-mem2',
        '--samtools', 'samtools',
        '--abra2', abra2_jar_path,
        '--abra2_targets_bed_file', abra2_targets_bed_file,
        '--gatk4', 'gatk',
        '--gatk4_baserecalibrator_params', '"--known-sites %s"' % known_sites_vcf_file,
        '--chromosomes', 'chr17',
        '--output_dir', output_dir,
        '-w', work_dir
    ]
    run_workflow(workflow='paired-end_read_dna_alignment_bwa-mem2.nf',
                 workflow_args=workflow_args)

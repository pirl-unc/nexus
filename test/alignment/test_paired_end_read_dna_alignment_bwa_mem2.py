import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path, get_alias_path


def test_paired_end_read_dna_alignment_bwa_mem2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    paired_end_read_tumor_dna_r1_fastq_file = get_data_path(name='fastq/hg38_tp53_tumor_paired-end_read_dna.r1.fastq.gz')
    paired_end_read_tumor_dna_r2_fastq_file = get_data_path(name='fastq/hg38_tp53_tumor_paired-end_read_dna.r2.fastq.gz')
    paired_end_read_normal_dna_r1_fastq_file = get_data_path(name='fastq/hg38_tp53_normal_paired-end_read_dna.r1.fastq.gz')
    paired_end_read_normal_dna_r2_fastq_file = get_data_path(name='fastq/hg38_tp53_normal_paired-end_read_dna.r2.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa')
    reference_genome_dict_file = get_data_path(name='fasta/hg38_chr17_1-8000000.dict')
    reference_genome_fasta_0123_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa.0123')
    reference_genome_fasta_amb_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa.amb')
    reference_genome_fasta_ann_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa.ann')
    reference_genome_fasta_bwt_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa.bwt.2bit.64')
    reference_genome_fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa.fai')
    reference_genome_fasta_pac_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa.pac')
    abra2_targets_bed_file = get_data_path(name='indices/abra2/gencode-v41-annotation-abra2-exon-targets.bed')
    known_sites_vcf_file = get_data_path(name='vcf/known_sites_hg38.vcf')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_paired_end_read_dna_alignment_bwa_mem2'
    work_dir = temp_dir + '/work/test_paired_end_read_dna_alignment_bwa_mem2'
    output_dir = temp_dir + '/outputs/test_paired_end_read_dna_alignment_bwa_mem2'
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
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_dict_file', reference_genome_dict_file,
        '--reference_genome_fasta_0123_file', reference_genome_fasta_0123_file,
        '--reference_genome_fasta_amb_file', reference_genome_fasta_amb_file,
        '--reference_genome_fasta_ann_file', reference_genome_fasta_ann_file,
        '--reference_genome_fasta_bwt_file', reference_genome_fasta_bwt_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--reference_genome_fasta_pac_file', reference_genome_fasta_pac_file,
        '--abra2_targets_bed_file', abra2_targets_bed_file,
        '--known_sites_files', known_sites_vcf_file,
        '--chromosomes', 'chr17',
        '--perform_local_indel_realignment', 'false',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='paired-end_read_dna_alignment_bwa-mem2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

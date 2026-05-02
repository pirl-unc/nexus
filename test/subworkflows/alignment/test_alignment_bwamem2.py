import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_alignment_bwamem2_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    paired_end_read_tumor_dna_r1_fastq_file = get_data_path(name='fastq/nexus-dna-001-tumor_paired-end_read_r1.fastq.gz')
    paired_end_read_tumor_dna_r2_fastq_file = get_data_path(name='fastq/nexus-dna-001-tumor_paired-end_read_r2.fastq.gz')
    paired_end_read_normal_dna_r1_fastq_file = get_data_path(name='fastq/nexus-dna-001-normal_paired-end_read_r1.fastq.gz')
    paired_end_read_normal_dna_r2_fastq_file = get_data_path(name='fastq/nexus-dna-001-normal_paired-end_read_r2.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    abra2_targets_bed_file = get_data_path(name='indices/abra2/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.abra2_targets.bed')
    known_sites_vcf_file = get_data_path(name='vcf/known_sites_hg38_chr17.vcf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_alignment_bwamem2'
    work_dir = temp_dir + '/work/test_alignment_bwamem2'
    output_dir = temp_dir + '/outputs/test_alignment_bwamem2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor', 'nexus-dna-001-normal'],
        'fastq_file_1': [paired_end_read_tumor_dna_r1_fastq_file, paired_end_read_normal_dna_r1_fastq_file],
        'fastq_file_2': [paired_end_read_tumor_dna_r2_fastq_file, paired_end_read_normal_dna_r2_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--abra2_targets_bed_file', abra2_targets_bed_file,
        '--known_sites_vcf_files', known_sites_vcf_file,
        '--chromosomes', '"chr17"',
        '--perform_local_indel_realignment', 'false',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='alignment_bwamem2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

def test_alignment_bwamem2_2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    paired_end_read_tumor_dna_r1_fastq_file = get_data_path(name='fastq/nexus-dna-001-tumor_paired-end_read_r1.fastq.gz')
    paired_end_read_tumor_dna_r2_fastq_file = get_data_path(name='fastq/nexus-dna-001-tumor_paired-end_read_r2.fastq.gz')
    paired_end_read_normal_dna_r1_fastq_file = get_data_path(name='fastq/nexus-dna-001-normal_paired-end_read_r1.fastq.gz')
    paired_end_read_normal_dna_r2_fastq_file = get_data_path(name='fastq/nexus-dna-001-normal_paired-end_read_r2.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    abra2_targets_bed_file = get_data_path(name='indices/abra2/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.abra2_targets.bed')
    known_sites_vcf_file = get_data_path(name='vcf/known_sites_hg38_chr17.vcf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_alignment_bwamem2_abra2'
    work_dir = temp_dir + '/work/test_alignment_bwamem2_abra2'
    output_dir = temp_dir + '/outputs/test_alignment_bwamem2_abra2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor', 'nexus-dna-001-normal'],
        'fastq_file_1': [paired_end_read_tumor_dna_r1_fastq_file, paired_end_read_normal_dna_r1_fastq_file],
        'fastq_file_2': [paired_end_read_tumor_dna_r2_fastq_file, paired_end_read_normal_dna_r2_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--abra2_targets_bed_file', abra2_targets_bed_file,
        '--known_sites_vcf_files', known_sites_vcf_file,
        '--abra2_temp_dir', os.getcwd() + '/tmp/abra2_temp/',
        '--chromosomes', '"chr17"',
        '--perform_local_indel_realignment', 'true',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='alignment_bwamem2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_alignment_bwamem2_multi_fastq_per_sample():
    """Two TSV rows share sample_id='nexus-dna-merged' but each row points to
    a different paired (R1, R2) fastq set. Both pairs should be concatenated
    in-process and aligned together into a single merged BAM via groupTuple
    — exercising the multi-row-per-sample-id code path. A second sample
    (single row) verifies the original one-pair-per-sample case still works
    alongside.
    """
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_r1 = get_data_path(name='fastq/nexus-dna-001-tumor_paired-end_read_r1.fastq.gz')
    tumor_r2 = get_data_path(name='fastq/nexus-dna-001-tumor_paired-end_read_r2.fastq.gz')
    normal_r1 = get_data_path(name='fastq/nexus-dna-001-normal_paired-end_read_r1.fastq.gz')
    normal_r2 = get_data_path(name='fastq/nexus-dna-001-normal_paired-end_read_r2.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    abra2_targets_bed_file = get_data_path(name='indices/abra2/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.abra2_targets.bed')
    known_sites_vcf_file = get_data_path(name='vcf/known_sites_hg38_chr17.vcf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_alignment_bwamem2_multi_fastq'
    work_dir = temp_dir + '/work/test_alignment_bwamem2_multi_fastq'
    output_dir = temp_dir + '/outputs/test_alignment_bwamem2_multi_fastq'
    for d in (intermediate_dir, work_dir, output_dir):
        if not os.path.exists(d):
            os.makedirs(d)
    # Two rows for 'nexus-dna-merged' (different paired fastq sets, same
    # sample_id), plus a single row for 'nexus-dna-001-normal' as a control.
    pd.DataFrame({
        'sample_id': [
            'nexus-dna-merged',
            'nexus-dna-merged',
            'nexus-dna-001-normal',
        ],
        'fastq_file_1': [tumor_r1, normal_r1, normal_r1],
        'fastq_file_2': [tumor_r2, normal_r2, normal_r2],
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--abra2_targets_bed_file', abra2_targets_bed_file,
        '--known_sites_vcf_files', known_sites_vcf_file,
        '--chromosomes', '"chr17"',
        '--perform_local_indel_realignment', 'false',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='alignment_bwamem2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

    # With perform_local_indel_realignment=false, the final BAM goes through
    # bwa-mem2 -> fixmate -> markdup -> ApplyBQSR, producing
    # ${sample_id}_bwamem2_sorted_fixmate_markeddup_recalibrated.bam.
    # Exactly one merged BAM should exist for the multi-row sample_id, and
    # one BAM for the single-row control.
    suffix = '_bwamem2_sorted_fixmate_markeddup_recalibrated.bam'
    merged_bam = os.path.join(output_dir, 'nexus-dna-merged' + suffix)
    merged_bai = merged_bam + '.bai'
    control_bam = os.path.join(output_dir, 'nexus-dna-001-normal' + suffix)
    control_bai = control_bam + '.bai'
    assert os.path.exists(merged_bam), \
        f"Expected merged BAM at {merged_bam}"
    assert os.path.exists(merged_bai), \
        f"Expected merged BAM index at {merged_bai}"
    assert os.path.exists(control_bam), \
        f"Expected control BAM at {control_bam}"
    assert os.path.exists(control_bai), \
        f"Expected control BAM index at {control_bai}"

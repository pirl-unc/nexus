import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path, get_alias_path


def test_paired_end_read_dna_variant_calling_svaba():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/sample100tumor_fixmate_markeddup_recalibrated.bam')
    tumor_bam_bai_file = get_data_path(name='bam/sample100tumor_fixmate_markeddup_recalibrated.bam.bai')
    normal_bam_file = get_data_path(name='bam/sample100normal_fixmate_markeddup_recalibrated.bam')
    normal_bam_bai_file = get_data_path(name='bam/sample100normal_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='indices/svaba/hg38_chr17_1-8M.fa')
    reference_genome_fasta_amb_file = get_data_path(name='indices/svaba/hg38_chr17_1-8M.fa.amb')
    reference_genome_fasta_ann_file = get_data_path(name='indices/svaba/hg38_chr17_1-8M.fa.ann')
    reference_genome_fasta_bwt_file = get_data_path(name='indices/svaba/hg38_chr17_1-8M.fa.bwt')
    reference_genome_fasta_fai_file = get_data_path(name='indices/svaba/hg38_chr17_1-8M.fa.fai')
    reference_genome_fasta_pac_file = get_data_path(name='indices/svaba/hg38_chr17_1-8M.fa.pac')
    reference_genome_fasta_sa_file = get_data_path(name='indices/svaba/hg38_chr17_1-8M.fa.sa')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_paired_end_read_dna_variant_calling_svaba'
    work_dir = temp_dir + '/work/test_paired_end_read_dna_variant_calling_svaba'
    output_dir = temp_dir + '/outputs/test_paired_end_read_dna_variant_calling_svaba'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample001tumor'],
        'tumor_bam_file': [tumor_bam_file],
        'tumor_bam_bai_file': [tumor_bam_bai_file],
        'normal_bam_file': [normal_bam_file],
        'normal_bam_bai_file': [normal_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_amb_file', reference_genome_fasta_amb_file,
        '--reference_genome_fasta_ann_file', reference_genome_fasta_ann_file,
        '--reference_genome_fasta_bwt_file', reference_genome_fasta_bwt_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--reference_genome_fasta_pac_file', reference_genome_fasta_pac_file,
        '--reference_genome_fasta_sa_file', reference_genome_fasta_sa_file,
        '--params_svaba', '"--hp --read-tracking"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='paired-end_read_dna_variant_calling_svaba.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


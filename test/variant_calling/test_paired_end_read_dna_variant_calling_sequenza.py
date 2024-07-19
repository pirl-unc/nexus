import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path, get_alias_path


def test_paired_end_read_dna_variant_calling_sequenza():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/hg38_chr21-22_tumor_paired_end_dna.bam')
    tumor_bam_bai_file = get_data_path(name='bam/hg38_chr21-22_tumor_paired_end_dna.bam.bai')
    normal_bam_file = get_data_path(name='bam/hg38_chr21-22_normal_paired_end_dna.bam')
    normal_bam_bai_file = get_data_path(name='bam/hg38_chr21-22_normal_paired_end_dna.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr21-22.fa.gz')
    reference_genome_wig_file = get_data_path(name='indices/sequenza/hg38_chr21-22.gc50.wig.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_paired_end_read_dna_variant_calling_sequenza'
    work_dir = temp_dir + '/work/test_paired_end_read_dna_variant_calling_sequenza'
    output_dir = temp_dir + '/outputs/test_paired_end_read_dna_variant_calling_sequenza'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample001'],
        'tumor_bam_file': [tumor_bam_file],
        'tumor_bam_bai_file': [tumor_bam_bai_file],
        'normal_bam_file': [normal_bam_file],
        'normal_bam_bai_file': [normal_bam_bai_file],
        'tumor_sample_id': ['tumor'],
        'normal_sample_id': ['normal']
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_wig_file', reference_genome_wig_file,
        '--chromosomes', '"chr21 chr22"',
        '--params_sequenza_bam2seqz', '"-N 5 --qformat sanger"',
        '--params_sequenza_seqzbinning', '"--window 50"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='paired-end_read_dna_variant_calling_sequenza.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


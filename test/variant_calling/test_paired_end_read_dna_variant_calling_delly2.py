import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path, get_alias_path


def test_paired_end_read_dna_variant_calling_delly2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/sample100tumor_fixmate_markeddup_recalibrated.bam')
    tumor_bam_bai_file = get_data_path(name='bam/sample100tumor_fixmate_markeddup_recalibrated.bam.bai')
    normal_bam_file = get_data_path(name='bam/sample100normal_fixmate_markeddup_recalibrated.bam')
    normal_bam_bai_file = get_data_path(name='bam/sample100normal_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa')
    reference_genome_fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.fai')
    exclude_tsv_file = get_data_path(name='indices/delly2/human.hg38.excl.tsv')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_paired_end_read_dna_variant_calling_delly2'
    work_dir = temp_dir + '/work/test_paired_end_read_dna_variant_calling_delly2'
    output_dir = temp_dir + '/outputs/test_paired_end_read_dna_variant_calling_delly2'
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
        'normal_bam_bai_file': [normal_bam_bai_file],
        'tumor_sample_id': ['tumor'],
        'normal_sample_id': ['normal']
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--exclude_tsv_file', exclude_tsv_file,
        '--params_delly2call', '""',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='paired-end_read_dna_variant_calling_delly2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


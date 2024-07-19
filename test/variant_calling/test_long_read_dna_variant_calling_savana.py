import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_long_read_dna_variant_calling_savana():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_dna_bam_file = get_data_path(name='bam/hg38_tp53_tumor_long_read_dna.bam')
    tumor_dna_bam_bai_file = get_data_path(name='bam/hg38_tp53_tumor_long_read_dna.bam.bai')
    normal_dna_bam_file = get_data_path(name='bam/hg38_tp53_normal_long_read_dna.bam')
    normal_dna_bam_bai_file = get_data_path(name='bam/hg38_tp53_normal_long_read_dna.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa')
    reference_genome_fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa.fai')
    savana_classification_configuration_file = get_data_path(name='indices/savana/savana_classification_parameters.json')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_dna_variant_calling_savana'
    work_dir = temp_dir + '/work/test_long_read_dna_variant_calling_savana'
    output_dir = temp_dir + '/outputs/test_long_read_dna_variant_calling_savana'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor'],
        'tumor_bam_file': [tumor_dna_bam_file],
        'tumor_bam_bai_file': [tumor_dna_bam_bai_file],
        'normal_bam_file': [normal_dna_bam_file],
        'normal_bam_bai_file': [normal_dna_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--custom_params_file', savana_classification_configuration_file,
        '--params_savana_run', '"--length 30 --mapq 20 --depth 3"',
        '--params_savana_classify', '""',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='long_read_dna_variant_calling_savana.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)
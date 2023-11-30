import pandas as pd
import os
from .data import get_data_path
from nexuslib.main import run_workflow


def test_paired_end_read_dna_somatic_small_variants_human():
    strelka2_path = os.environ.get('strelka2')
    python2_path = os.environ.get('python2')
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    tumor_bam_file = get_data_path(name='hg38_tp53_tumor_paired_end_dna_fixmate_markeddup_recalibrated.bam')
    tumor_bam_bai_file = get_data_path(name='hg38_tp53_tumor_paired_end_dna_fixmate_markeddup_recalibrated.bam.bai')
    normal_bam_file = get_data_path(name='hg38_tp53_normal_paired_end_dna_fixmate_markeddup_recalibrated.bam')
    normal_bam_bai_file = get_data_path(name='hg38_tp53_normal_paired_end_dna_fixmate_markeddup_recalibrated.bam.bai')
    known_sites_vcf_file = get_data_path(name='known_sites_hg38.vcf')
    reference_genome_fasta_file = get_data_path(name='hg38_chr17_1-8000000.fa')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_paired_end_read_dna_somatic_small_variants_human'
    work_dir = temp_dir + '/work/test_paired_end_read_dna_somatic_small_variants_human'
    output_dir = temp_dir + '/outputs/test_paired_end_read_dna_somatic_small_variants_human'
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
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--is_human', 'true',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--python2', python2_path,
        '--gatk4', 'gatk',
        '--gatk4_mutect2_params', '"--germline-resource %s --panel-of-normals %s "' % (known_sites_vcf_file, known_sites_vcf_file),
        '--gatk4_getpileupsummaries_params', '"-V %s -L %s "' % (known_sites_vcf_file, known_sites_vcf_file),
        '--gatk4_chromosomes', 'chr17',
        '--picard', 'picard',
        '--strelka2', strelka2_path,
        '--strelka2_params', " ",
        '--containerization', "docker",
        '--deepvariant_bin_path', '/opt/deepvariant/bin/run_deepvariant',
        '--deepvariant_input_path', '/Users/leework/Documents/Research/projects/nexus/',
        '--deepvariant_output_path', '/var/folders/',
        '--deepvariant_bin_version', '1.6.0',
        '--deepvariant_model_type', 'WGS',
        '--output_dir', output_dir,
        '-w', work_dir
    ]
    run_workflow(workflow='paired-end_read_dna_somatic_small_variants.nf',
                 workflow_args=workflow_args)

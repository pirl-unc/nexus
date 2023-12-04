import pandas as pd
import os
from .data import get_data_path, get_alias_path
from nexuslib.main import run_workflow


def test_paired_end_read_dna_somatic_structural_variants():
    delly2_path = get_alias_path(executable_name='delly')
    python2_path = get_alias_path(executable_name='python2')
    lumpyexpress_path = get_alias_path(executable_name='lumpyexpress')
    lumpyexpress_config_path = get_alias_path(executable_name='lumpyexpress.config')
    lumpy_extract_split_reads_script_file_path = get_alias_path(executable_name='extractSplitReads_BwaMem')
    samtools_path = get_alias_path(executable_name='samtools')
    bcftools_path = get_alias_path(executable_name='bcftools')
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    tumor_bam_file = get_data_path(name='hg38_tp53_tumor_paired_end_dna_fixmate_markeddup_recalibrated.bam')
    tumor_bam_bai_file = get_data_path(name='hg38_tp53_tumor_paired_end_dna_fixmate_markeddup_recalibrated.bam.bai')
    normal_bam_file = get_data_path(name='hg38_tp53_normal_paired_end_dna_fixmate_markeddup_recalibrated.bam')
    normal_bam_bai_file = get_data_path(name='hg38_tp53_normal_paired_end_dna_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='hg38_chr17_1-8000000.fa')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_paired_end_read_dna_somatic_structural_variants'
    work_dir = temp_dir + '/work/test_paired_end_read_dna_somatic_structural_variants'
    output_dir = temp_dir + '/outputs/test_paired_end_read_dna_somatic_structural_variants'
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
        '--delly2', delly2_path,
        '--delly2_params', '"--map-qual 20 "',
        '--bcftools', bcftools_path,
        '--python2', python2_path,
        '--lumpyexpress', lumpyexpress_path,
        '--lumpyexpress_config_file', lumpyexpress_config_path,
        '--lumpy_extract_split_reads_script_file', lumpy_extract_split_reads_script_file_path,
        '--samtools', samtools_path,
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='paired-end_read_dna_somatic_structural_variants.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

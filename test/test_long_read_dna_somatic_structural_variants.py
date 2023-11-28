# import pandas as pd
# import os
# from .data import get_data_path
# from nexuslib.main import run_workflow
#
#
# def test_pacbio_dna_somatic_structural_variants():
#     nextflow_config_file = get_data_path(name='nextflow_test.config')
#     long_read_tumor_dna_bam_file = get_data_path(name='hg38_tp53_tumor_dna_minimap2_mdtagged_sorted.bam')
#     long_read_tumor_dna_bam_bai_file = get_data_path(name='hg38_tp53_tumor_dna_minimap2_mdtagged_sorted.bam.bai')
#     long_read_normal_dna_bam_file = get_data_path(name='hg38_tp53_normal_dna_minimap2_mdtagged_sorted.bam')
#     long_read_normal_dna_bam_bai_file = get_data_path(name='hg38_tp53_normal_dna_minimap2_mdtagged_sorted.bam.bai')
#     reference_genome_fasta_file = get_data_path(name='hg38_chr17_1-8000000.fa')
#     temp_dir = os.getcwd() + '/tmp/'
#     intermediate_dir = temp_dir + 'intermediate/test_pacbio_dna_somatic_structural_variants/'
#     work_dir = temp_dir + 'work/test_pacbio_dna_somatic_structural_variants/'
#     output_dir = temp_dir + 'outputs/test_pacbio_dna_somatic_structural_variants/'
#     if not os.path.exists(intermediate_dir):
#         os.makedirs(intermediate_dir)
#     if not os.path.exists(work_dir):
#         os.makedirs(work_dir)
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
#     pd.DataFrame({
#         'sample_id': ['sample001'],
#         'tumor_bam_file': [long_read_tumor_dna_bam_file],
#         'tumor_bam_bai_file': [long_read_tumor_dna_bam_bai_file],
#         'normal_bam_file': [long_read_normal_dna_bam_file],
#         'normal_bam_bai_file': [long_read_normal_dna_bam_bai_file]
#     }).to_csv(intermediate_dir + "samples.tsv", sep='\t', index=False)
#     workflow_args = [
#         '-c', nextflow_config_file,
#         '--samples_tsv_file', intermediate_dir + 'samples.tsv',
#         '--reference_genome_fasta_file', reference_genome_fasta_file,
#         '--savana', 'savana',
#         '--savana_params', '"--length 30 --mapq 20 --depth 3 "',
#         '--output_dir', output_dir,
#         '-w', work_dir
#     ]
#     run_workflow(workflow='pacbio_dna_somatic_structural_variants',
#                  workflow_args=workflow_args)

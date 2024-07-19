import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_long_read_dna_variant_calling_cutesv():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_dna_bam_file = get_data_path(name='bam/hg38_tp53_tumor_long_read_dna.bam')
    tumor_dna_bam_bai_file = get_data_path(name='bam/hg38_tp53_tumor_long_read_dna.bam.bai')
    normal_dna_bam_file = get_data_path(name='bam/hg38_tp53_normal_long_read_dna.bam')
    normal_dna_bam_bai_file = get_data_path(name='bam/hg38_tp53_normal_long_read_dna.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_dna_variant_calling_cutesv'
    work_dir = temp_dir + '/work/test_long_read_dna_variant_calling_cutesv'
    output_dir = temp_dir + '/outputs/test_long_read_dna_variant_calling_cutesv'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor', 'normal'],
        'bam_file': [tumor_dna_bam_file, normal_dna_bam_file],
        'bam_bai_file': [tumor_dna_bam_bai_file, normal_dna_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--params_cutesv', '"--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_support 3 --min_mapq 20 --min_size 30 --max_size -1 --report_readid --genotype"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='long_read_dna_variant_calling_cutesv.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

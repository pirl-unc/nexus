import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_long_read_dna_variant_phasing_whatshap():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_dna_bam_file = get_data_path(name='bam/sample001tumor_minimap2_mdtagged_sorted.bam')
    long_read_tumor_dna_bam_bai_file = get_data_path(name='bam/sample001tumor_minimap2_mdtagged_sorted.bam.bai')
    small_variants_vcf_gz_file = get_data_path(name='vcf/sample001tumor_deepvariant.vcf')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa')
    reference_genome_fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa.fai')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_dna_variant_phasing_whatshap'
    work_dir = temp_dir + '/work/test_long_read_dna_variant_phasing_whatshap'
    output_dir = temp_dir + '/outputs/test_long_read_dna_variant_phasing_whatshap'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample001tumor'],
        'bam_file': [long_read_tumor_dna_bam_file],
        'bam_bai_file': [long_read_tumor_dna_bam_bai_file],
        'small_variants_vcf_file': [small_variants_vcf_gz_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='long_read_dna_variant_phasing_whatshap.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

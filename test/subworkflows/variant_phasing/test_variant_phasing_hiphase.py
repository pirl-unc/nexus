import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_phasing_hiphase():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_dna_bam_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    long_read_tumor_dna_bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    small_variants_vcf_gz_file = get_data_path(name='vcf/nexus-dna-002-tumor_whatshap_phased.vcf.gz')
    small_variants_vcf_gz_tbi_file = get_data_path(name='vcf/nexus-dna-002-tumor_whatshap_phased.vcf.gz.tbi')
    structural_variants_vcf_gz_file =get_data_path(name='vcf/nexus-dna-002-tumor_pbsv.vcf.gz')
    structural_variants_vcf_gz_tbi_file =get_data_path(name='vcf/nexus-dna-002-tumor_pbsv.vcf.gz.tbi')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_phasing_hiphase'
    work_dir = temp_dir + '/work/test_variant_phasing_hiphase'
    output_dir = temp_dir + '/outputs/test_variant_phasing_hiphase'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-002-tumor-long-read'],
        'bam_file': [long_read_tumor_dna_bam_file],
        'bam_bai_file': [long_read_tumor_dna_bam_bai_file],
        'small_variants_vcf_gz_file': [small_variants_vcf_gz_file],
        'small_variants_vcf_gz_tbi_file': [small_variants_vcf_gz_tbi_file],
        'structural_variants_vcf_gz_file': [structural_variants_vcf_gz_file],
        'structural_variants_vcf_gz_tbi_file': [structural_variants_vcf_gz_tbi_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_phasing_hiphase.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

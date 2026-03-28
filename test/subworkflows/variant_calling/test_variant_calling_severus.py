import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_calling_severus():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_dna_bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    tumor_dna_bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    normal_dna_bam_file = get_data_path(name='bam/nexus-dna-001-normal-long-read_minimap2_mdtagged_sorted.bam')
    normal_dna_bam_bai_file = get_data_path(name='bam/nexus-dna-001-normal-long-read_minimap2_mdtagged_sorted.bam.bai')
    phased_vcf_file = get_data_path(name='vcf/nexus-dna-001-tumor_whatshap_phased.vcf.gz')
    vntr_bed_file = get_data_path(name='indices/severus/human_GRCh38_no_alt_analysis_set.trf.bed')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_severus'
    work_dir = temp_dir + '/work/test_variant_calling_severus'
    output_dir = temp_dir + '/outputs/test_variant_calling_severus'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-002-tumor'],
        'tumor_bam_file': [tumor_dna_bam_file],
        'tumor_bam_bai_file': [tumor_dna_bam_bai_file],
        'normal_bam_file': [normal_dna_bam_file],
        'normal_bam_bai_file': [normal_dna_bam_bai_file],
        'phased_vcf_file': [phased_vcf_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--vntr_bed_file', vntr_bed_file,
        '--params_severus', '"--min-support 3 --min-sv-size 10 --min-mapq 20 --output-read-ids --bp-cluster-size 50"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_severus.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_haplotagging_margin():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_dna_bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    long_read_tumor_dna_bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    small_variants_vcf_gz_file = get_data_path(name='vcf/nexus-dna-001-tumor_deepvariant.vcf.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    margin_phase_params_json_file = get_data_path(
        name='indices/margin/phase/allParams.phase_vcf.pb-hifi.json'
    )
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_margin'
    work_dir = temp_dir + '/work/test_haplotagging_margin'
    output_dir = temp_dir + '/outputs/test_haplotagging_margin'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'bam_file': [long_read_tumor_dna_bam_file],
        'bam_bai_file': [long_read_tumor_dna_bam_bai_file],
        'small_variants_vcf_file': [small_variants_vcf_gz_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--output_dir', output_dir,
        '--margin_phase_params_json_file', margin_phase_params_json_file,
    ]
    run_workflow(workflow='haplotagging_margin.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

    # Verify both outputs were produced.
    sample_out_dir = os.path.join(output_dir, 'nexus-dna-001-tumor')
    phased_vcf = os.path.join(sample_out_dir, 'nexus-dna-001-tumor_margin_phased.vcf.gz')
    phased_vcf_tbi = os.path.join(sample_out_dir, 'nexus-dna-001-tumor_margin_phased.vcf.gz.tbi')
    haplotagged_bam = os.path.join(
        sample_out_dir,
        'nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted_haplotagged.bam'
    )
    haplotagged_bai = haplotagged_bam + '.bai'
    assert os.path.exists(phased_vcf), f"Expected phased VCF at {phased_vcf}"
    assert os.path.exists(phased_vcf_tbi), f"Expected phased VCF index at {phased_vcf_tbi}"
    assert os.path.exists(haplotagged_bam), f"Expected haplotagged BAM at {haplotagged_bam}"
    assert os.path.exists(haplotagged_bai), f"Expected haplotagged BAM index at {haplotagged_bai}"

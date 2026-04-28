import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_haplotype_phasing_hapcut2_whatshap():
    """End-to-end test for the chained HapCUT2 + WhatsHap haplotag subworkflow.

    Runs `extractHAIRS` + `HAPCUT2` to produce a phased VCF, then runs
    `whatshap haplotag` to produce a haplotagged BAM.
    """
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_dna_bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    long_read_tumor_dna_bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    small_variants_vcf_gz_file = get_data_path(name='vcf/nexus-dna-001-tumor_deepvariant.vcf.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotype_phasing_hapcut2_whatshap'
    work_dir = temp_dir + '/work/test_haplotype_phasing_hapcut2_whatshap'
    output_dir = temp_dir + '/outputs/test_haplotype_phasing_hapcut2_whatshap'
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
        '--read_technology', 'pacbio',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='haplotype_phasing_hapcut2-whatshap.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_calling_savana():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_dna_bam_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    tumor_dna_bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    normal_dna_bam_file = get_data_path(name='bam/nexus-dna-002-normal-long-read_minimap2_mdtagged_sorted.bam')
    normal_dna_bam_bai_file = get_data_path(name='bam/nexus-dna-002-normal-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    savana_classification_configuration_file = get_data_path(name='indices/savana/savana_classification_parameters.json')
    contigs_file = get_data_path(name='indices/savana/hg38_chr17_contig.txt')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_savana'
    work_dir = temp_dir + '/work/test_variant_calling_savana'
    output_dir = temp_dir + '/outputs/test_variant_calling_savana'
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
        'normal_bam_bai_file': [normal_dna_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--contigs_txt_file', contigs_file,
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--custom_params_file', savana_classification_configuration_file,
        '--params_savana_run', '"--length 10 --mapq 20 --min_support 3"',
        '--params_savana_classify', '""',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_savana.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

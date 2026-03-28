import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_calling_reditools2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    rna_bam_file = get_data_path(name='bam/nexus-rna-101-tumor-paired-end-read_star_outputs/nexus-rna-101-tumor-paired-end-read_star_Aligned.sortedByCoord.out.bam')
    rna_bam_bai_file = get_data_path(name='bam/nexus-rna-101-tumor-paired-end-read_star_outputs/nexus-rna-101-tumor-paired-end-read_star_Aligned.sortedByCoord.out.bam.bai')
    dna_bam_file = get_data_path(name='bam/nexus-dna-101-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    dna_bam_bai_file = get_data_path(name='bam/nexus-dna-101-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.sorted.filtered.gtf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_reditools2'
    work_dir = temp_dir + '/work/test_variant_calling_reditools2'
    output_dir = temp_dir + '/outputs/test_variant_calling_reditools2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-101-tumor'],
        'rna_bam_file': [rna_bam_file],
        'rna_bam_bai_file': [rna_bam_bai_file],
        'dna_bam_file': [dna_bam_file],
        'dna_bam_bai_file': [dna_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--params_reditools_annotatetable', '"-s 4 -c 1,2,3 -n gencode"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_reditools2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

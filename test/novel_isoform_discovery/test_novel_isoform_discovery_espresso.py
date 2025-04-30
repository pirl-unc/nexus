import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_novel_isoform_discovery_espresso():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_bam_file = get_data_path(name='bam/sample200tumor_minimap2_mdtagged_sorted.bam')
    long_read_tumor_bam_bai_file = get_data_path(name='bam/sample200tumor_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa')
    reference_genome_fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa.fai')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode_v41_tp53_annotation.gtf')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_novel_isoform_discovery_espresso'
    work_dir = temp_dir + '/work/test_novel_isoform_discovery_espresso'
    output_dir = temp_dir + '/outputs/test_novel_isoform_discovery_espresso'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample200tumor'],
        'bam_file': [long_read_tumor_bam_file],
        'bam_bai_file': [long_read_tumor_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--output_dir', output_dir
    ]
    run_workflow(workflow='novel_isoform_discovery_espresso.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

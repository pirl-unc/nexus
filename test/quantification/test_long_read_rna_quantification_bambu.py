import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path, get_alias_path


def test_long_read_rna_quantification_bambu():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/sample200tumor_minimap2_mdtagged_sorted.bam')
    tumor_bam_bai_file = get_data_path(name='bam/sample200tumor_minimap2_mdtagged_sorted.bam.bai')
    normal_bam_file = get_data_path(name='bam/sample200normal_minimap2_mdtagged_sorted.bam')
    normal_bam_bai_file = get_data_path(name='bam/sample200normal_minimap2_mdtagged_sorted.bam.bai')
    fasta_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa')
    fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8M.fa.fai')
    gtf_file = get_data_path(name='gtf/gencode_v41_tp53_annotation.gtf')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_rna_quantification_bambu'
    work_dir = temp_dir + '/work/test_long_read_rna_quantification_bambu'
    output_dir = temp_dir + '/outputs/test_long_read_rna_quantification_bambu'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample200tumor', 'sample200normal'],
        'bam_file': [tumor_bam_file, normal_bam_file],
        'bam_bai_file': [tumor_bam_bai_file, normal_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', fasta_file,
        '--reference_genome_fasta_fai_file', fasta_fai_file,
        '--gtf_file', gtf_file,
        '--output_dir', output_dir
    ]
    run_workflow(workflow='long_read_rna_quantification_bambu.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

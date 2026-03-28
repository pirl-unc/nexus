import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_isoform_characterization_sqanti3_gtf():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_tumor_rna_gtf_file = get_data_path(name='inputs/subworkflows/sqanti3/nexus-rna-002-tumor_flair.isoforms.gtf')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_isoform_characterization_sqanti3_gtf'
    work_dir = temp_dir + '/work/test_isoform_characterization_sqanti3_gtf'
    output_dir = temp_dir + '/outputs/test_isoform_characterization_sqanti3_gtf'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-002-tumor_flair'],
        'gtf_file': [long_read_tumor_rna_gtf_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--params_sqanti3_qc', '"--report skip"',
        '--params_sqanti3_filter', '"ml --skip_report"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='isoform_characterization_sqanti3-gtf.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

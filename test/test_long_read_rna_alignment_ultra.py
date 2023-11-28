import pandas as pd
import os
from .data import get_data_path
from nexuslib.main import run_workflow


def test_long_read_rna_alignment_ultra():
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    long_read_tumor_rna_fastq_file = get_data_path(name='hg38_tp53_tumor_rna.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='hg38_chr17_1-8000000.fa')
    reference_gtf_file = get_data_path(name='gencode_v41_tp53_annotation.gtf')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_rna_alignment_ultra'
    work_dir = temp_dir + '/work/test_long_read_rna_alignment_ultra/'
    output_dir = temp_dir + '/outputs/test_long_read_rna_alignment_ultra/'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor'],
        'fastq_file': [long_read_tumor_rna_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    # Create index
    if not os.path.exists(intermediate_dir + '/uLTRA_index/'):
        os.makedirs(intermediate_dir + '/uLTRA_index/')
    cmd = [
        'uLTRA', 'index',
        '--disable_infer',
        reference_genome_fasta_file,
        reference_gtf_file,
        intermediate_dir + '/uLTRA_index/'
    ]
    os.system(' '.join(cmd))

    # Run uLTRA
    workflow_args = [
        '-c', nextflow_config_file,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--ultra', 'uLTRA',
        '--ultra_index', intermediate_dir + '/uLTRA_index/',
        '--ultra_params', '"--isoseq "', # keep the space at the end
        '--samtools', 'samtools',
        '--output_dir', output_dir,
        '-w', work_dir
    ]
    run_workflow(workflow='long_read_rna_alignment_ultra.nf',
                 workflow_args=workflow_args)

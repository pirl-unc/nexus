import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_paired_end_read_rna_alignment_star():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    paired_end_read_tumor_rna_r1_fastq_file = get_data_path(name='fastq/sample301tumor_paired-end_read_rna.r1.fastq.gz')
    paired_end_read_tumor_rna_r2_fastq_file = get_data_path(name='fastq/sample301tumor_paired-end_read_rna.r2.fastq.gz')
    star_index = get_data_path(name='indices/star/hg38_chr17/')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_paired_end_read_rna_alignment_star'
    work_dir = temp_dir + '/work/test_paired_end_read_rna_alignment_star'
    output_dir = temp_dir + '/outputs/test_paired_end_read_rna_alignment_star'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample301tumor'],
        'fastq_file_1': [paired_end_read_tumor_rna_r1_fastq_file],
        'fastq_file_2': [paired_end_read_tumor_rna_r2_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--star_index', star_index,
        '--params_star', '"--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='paired-end_read_rna_alignment_star.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


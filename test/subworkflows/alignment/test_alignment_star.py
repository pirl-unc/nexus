import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_alignment_star():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_1_fastq_file_1 = get_data_path(name='fastq/nexus-rna-001-tumor_paired-end_read_r1.fastq.gz')
    tumor_1_fastq_file_2 = get_data_path(name='fastq/nexus-rna-001-tumor_paired-end_read_r2.fastq.gz')
    tumor_2_fastq_file_1 = get_data_path(name='fastq/nexus-rna-004-tumor_paired-end_read_r1.fastq.gz')
    tumor_2_fastq_file_2 = get_data_path(name='fastq/nexus-rna-004-tumor_paired-end_read_r2.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_alignment_star'
    work_dir = temp_dir + '/work/test_alignment_star'
    output_dir = temp_dir + '/outputs/test_alignment_star'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-001-tumor', 'nexus-rna-004-tumor'],
        'fastq_file_1': [tumor_1_fastq_file_1, tumor_2_fastq_file_1],
        'fastq_file_2': [tumor_1_fastq_file_2, tumor_2_fastq_file_2]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--params_star_genomegenerate', '"--genomeSAindexNbases 10 --sjdbOverhang 150"',
        '--params_star', '"--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic --chimSegmentMin 10"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='alignment_star.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path, get_alias_path


def test_paired_end_read_rna_variant_calling_arriba():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fastq_file_1 = get_data_path(name='fastq/hg38_tp53_tumor_paired-end_read_rna.r1.fastq.gz')
    fastq_file_2 = get_data_path(name='fastq/hg38_tp53_tumor_paired-end_read_rna.r2.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa')
    reference_genome_fasta_fai_file = get_data_path(name='fasta/hg38_chr17_1-8000000.fa.fai')
    star_index = get_data_path(name='indices/star/hg38_chr17/')
    gtf_file = get_data_path(name='gtf/gencode_v41_tp53_annotation.gtf')
    protein_domains_gff3_file = get_data_path(name='indices/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_paired_end_read_rna_variant_calling_arriba'
    work_dir = temp_dir + '/work/test_paired_end_read_rna_variant_calling_arriba'
    output_dir = temp_dir + '/outputs/test_paired_end_read_rna_variant_calling_arriba'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample001'],
        'fastq_file_1': [fastq_file_1],
        'fastq_file_2': [fastq_file_2]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genome_fasta_fai_file', reference_genome_fasta_fai_file,
        '--star_index', star_index,
        '--gtf_file', gtf_file,
        '--protein_domains_gff3_file', protein_domains_gff3_file,
        '--params_arriba', '"-S 3 -f blacklist -i chr*"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='paired-end_read_rna_variant_calling_arriba.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


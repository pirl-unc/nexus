import pandas as pd
import os
from .data import get_data_path
from nexuslib.main import run_workflow


def test_novel_isoform_discovery_isoquant():
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    long_read_tumor_rna_fastq_file = get_data_path(name='hg38_tp53_tumor_long_read_rna.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='hg38_chr17_1-8000000.fa')
    gtf_file = get_data_path(name='gencode_v41_tp53_annotation.gtf')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_novel_isoform_discovery_isoquant'
    work_dir = temp_dir + '/work/test_novel_isoform_discovery_isoquant'
    output_dir = temp_dir + '/outputs/test_novel_isoform_discovery_isoquant'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(intermediate_dir + '/fastq_list.txt', 'w') as f:
        f.write('#EXPERIMENT1\n')
        f.write('%s:%s\n' % (long_read_tumor_rna_fastq_file, 'tumor'))
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--isoquant_params', '"--data_type pacbio_ccs --complete_genedb --high_memory --genedb %s --fastq_list %s "' % (gtf_file, intermediate_dir + '/fastq_list.txt'),
        '--output_dir', output_dir
    ]
    run_workflow(workflow='novel_isoform_discovery_isoquant.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


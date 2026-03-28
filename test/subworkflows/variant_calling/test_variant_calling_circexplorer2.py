import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_calling_circexplorer2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    input_file = get_data_path(name='bam/nexus-rna-102-tumor-paired-end-read_star_outputs/nexus-rna-102-tumor-paired-end-read_star_Chimeric.out.junction')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    circexplorer2_gene_annotation_txt_file = get_data_path(name='indices/circexplorer2/hg38_kg/hg38_kg.txt')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_circexplorer2'
    work_dir = temp_dir + '/work/test_variant_calling_circexplorer2'
    output_dir = temp_dir + '/outputs/test_variant_calling_circexplorer2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-102-tumor'],
        'input_file': [input_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--circexplorer2_gene_annotation_txt_file', circexplorer2_gene_annotation_txt_file,
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_circexplorer2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_calling_arriba():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-rna-004-tumor-paired-end-read_star_outputs/nexus-rna-004-tumor-paired-end-read_star_Aligned.sortedByCoord.out.bam')
    bai_file = get_data_path(name='bam/nexus-rna-004-tumor-paired-end-read_star_outputs/nexus-rna-004-tumor-paired-end-read_star_Aligned.sortedByCoord.out.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    protein_domains_gff3_file = get_data_path(name='indices/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_arriba'
    work_dir = temp_dir + '/work/test_variant_calling_arriba'
    output_dir = temp_dir + '/outputs/test_variant_calling_arriba'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-004-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--reference_genes_gtf_file', reference_genes_gtf_file,
        '--protein_domains_gff3_file', protein_domains_gff3_file,
        '--params_arriba', '"-S 3 -f blacklist -i chr*"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_arriba.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

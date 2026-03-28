import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_peptide_prediction_mopepgen():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    mutect2_vcf_file = get_data_path(name='vcf/nexus-dna-001-tumor_gatk4-mutect2.vcf')
    strelka2_snv_vcf_file = get_data_path(name='vcf/nexus-dna-001-tumor_strelka2_snvs.vcf')
    strelka2_indel_vcf_file = get_data_path(name='vcf/nexus-dna-002-tumor_strelka2_indels.vcf')
    reditools2_tsv_file = get_data_path(name='tsv/nexus-rna-101-tumor_rna_reditools2_final_annotated.tsv')
    arriba_tsv_file = get_data_path(name='tsv/nexus-rna-004-tumor_arriba.tsv')
    rmats_output_dir = get_data_path(name='inputs/subworkflows/mopepgen/nexus-rna-103-tumor_rmats_outputs/')
    circexplorer2_tsv_file = get_data_path(name='tsv/nexus-rna-102-tumor_circexplorer2_back_spliced_junction_annotated.bed')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_genes_gtf_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    reference_proteome_fasta_file = get_data_path(name='fasta/gencode.v45.pc_translations.fa')
    params_yaml_file = get_data_path(name='inputs/subworkflows/peptide_prediction_mopepgen-gencode/params.yaml')

    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_peptide_prediction_mopepgen'
    work_dir = temp_dir + '/work/test_peptide_prediction_mopepgen'
    output_dir = temp_dir + '/outputs/test_peptide_prediction_mopepgen'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample1'],
        'mutect2_vcf_file': [mutect2_vcf_file],
        'strelka2_snv_vcf_file': [strelka2_snv_vcf_file],
        'strelka2_indel_vcf_file': [strelka2_indel_vcf_file],
        'reditools2_tsv_file': [reditools2_tsv_file],
        'arriba_tsv_file': [arriba_tsv_file],
        'rmats_output_dir': [rmats_output_dir],
        'circexplorer2_tsv_file': [circexplorer2_tsv_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['reference_genes_gtf_file'] = reference_genes_gtf_file
    params['reference_proteome_fasta_file'] = reference_proteome_fasta_file
    params['vep']['dir'] = '/Users/leework/.vep'
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file
    ]
    run_workflow(workflow='peptide_prediction_mopepgen-gencode.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

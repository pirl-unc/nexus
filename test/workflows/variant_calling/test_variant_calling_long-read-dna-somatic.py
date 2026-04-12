import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


# ClairS, DeepSomatic (SNV/indel callers — nexus-dna-001 tumor/normal)
def test_variant_calling_long_read_dna_somatic_github_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    tumor_bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    normal_bam_file = get_data_path(name='bam/nexus-dna-001-normal-long-read_minimap2_mdtagged_sorted.bam')
    normal_bam_bai_file = get_data_path(name='bam/nexus-dna-001-normal-long-read_minimap2_mdtagged_sorted.bam.bai')
    normal_small_variants_vcf_file = get_data_path(name='vcf/nexus-dna-001-tumor_deepvariant.vcf.gz')  # using tumor DV VCF as stand-in for normal
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/variant_calling_long-read-dna-somatic/params_1.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_longread_somatic_1'
    work_dir = temp_dir + '/work/test_variant_calling_longread_somatic_1'
    output_dir = temp_dir + '/outputs/test_variant_calling_longread_somatic_1'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'tumor_bam_file': [tumor_bam_file],
        'tumor_bam_bai_file': [tumor_bam_bai_file],
        'normal_bam_file': [normal_bam_file],
        'normal_bam_bai_file': [normal_bam_bai_file],
        'normal_small_variants_vcf_file': [normal_small_variants_vcf_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['deepsomatic']['input_path'] = '/home/runner/work/nexus/nexus/'
    params['deepsomatic']['output_path'] = '/tmp/'
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='variant_calling_long-read-dna-somatic.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# ClairS, DeepSomatic (SNV/indel callers — nexus-dna-001 tumor/normal)
def test_variant_calling_long_read_dna_somatic_local_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    tumor_bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    normal_bam_file = get_data_path(name='bam/nexus-dna-001-normal-long-read_minimap2_mdtagged_sorted.bam')
    normal_bam_bai_file = get_data_path(name='bam/nexus-dna-001-normal-long-read_minimap2_mdtagged_sorted.bam.bai')
    normal_small_variants_vcf_file = get_data_path(name='vcf/nexus-dna-001-tumor_deepvariant.vcf.gz')  # using tumor DV VCF as stand-in for normal
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/variant_calling_long-read-dna-somatic/params_1.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_longread_somatic_1'
    work_dir = temp_dir + '/work/test_variant_calling_longread_somatic_1'
    output_dir = temp_dir + '/outputs/test_variant_calling_longread_somatic_1'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'tumor_bam_file': [tumor_bam_file],
        'tumor_bam_bai_file': [tumor_bam_bai_file],
        'normal_bam_file': [normal_bam_file],
        'normal_bam_bai_file': [normal_bam_bai_file],
        'normal_small_variants_vcf_file': [normal_small_variants_vcf_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['deepsomatic']['input_path'] = '/Users/leework/Documents/Research/projects/project_nexus/'
    params['deepsomatic']['output_path'] = '/var/folders/'
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='variant_calling_long-read-dna-somatic.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# Delly2, Nanomonsv (SV callers — nexus-dna-005 tumor/normal)
def test_variant_calling_long_read_dna_somatic_2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/nexus-dna-005-tumor-long-read_minimap2_mdtagged_sorted.bam')
    tumor_bam_bai_file = get_data_path(name='bam/nexus-dna-005-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    normal_bam_file = get_data_path(name='bam/nexus-dna-005-normal-long-read_minimap2_mdtagged_sorted.bam')
    normal_bam_bai_file = get_data_path(name='bam/nexus-dna-005-normal-long-read_minimap2_mdtagged_sorted.bam.bai')
    normal_small_variants_vcf_file = get_data_path(name='vcf/nexus-dna-005-tumor_deepvariant.vcf.gz')  # placeholder; Severus not run in this test
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr21-22.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/variant_calling_long-read-dna-somatic/params_2.yaml')
    delly2_exclude_tsv_file = get_data_path(name='indices/delly2/human.hg38.excl.tsv')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_longread_somatic_2'
    work_dir = temp_dir + '/work/test_variant_calling_longread_somatic_2'
    output_dir = temp_dir + '/outputs/test_variant_calling_longread_somatic_2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-005-tumor'],
        'tumor_bam_file': [tumor_bam_file],
        'tumor_bam_bai_file': [tumor_bam_bai_file],
        'normal_bam_file': [normal_bam_file],
        'normal_bam_bai_file': [normal_bam_bai_file],
        'normal_small_variants_vcf_file': [normal_small_variants_vcf_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['delly2']['exclude_tsv_file'] = delly2_exclude_tsv_file

    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='variant_calling_long-read-dna-somatic.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# Savana, Severus, SVision-pro (SV callers — nexus-dna-002 tumor/normal)
def test_variant_calling_long_read_dna_somatic_3():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    tumor_bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    normal_bam_file = get_data_path(name='bam/nexus-dna-002-normal-long-read_minimap2_mdtagged_sorted.bam')
    normal_bam_bai_file = get_data_path(name='bam/nexus-dna-002-normal-long-read_minimap2_mdtagged_sorted.bam.bai')
    normal_small_variants_vcf_file = get_data_path(name='vcf/nexus-dna-002-tumor_deepvariant.vcf.gz')  # using tumor DV VCF as stand-in for normal
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/variant_calling_long-read-dna-somatic/params_3.yaml')
    savana_contigs_txt_file = get_data_path(name='indices/savana/hg38_chr17_contig.txt')
    savana_custom_params_file = get_data_path(name='indices/savana/savana_classification_parameters.json')
    severus_vntr_bed_file = get_data_path(name='indices/severus/human_GRCh38_no_alt_analysis_set.trf.bed')
    svisionpro_model_file = get_data_path(name='indices/svisionpro/model_liteunet_256_8_16_32_32_32.pth')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_longread_somatic_3'
    work_dir = temp_dir + '/work/test_variant_calling_longread_somatic_3'
    output_dir = temp_dir + '/outputs/test_variant_calling_longread_somatic_3'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-002-tumor'],
        'tumor_bam_file': [tumor_bam_file],
        'tumor_bam_bai_file': [tumor_bam_bai_file],
        'normal_bam_file': [normal_bam_file],
        'normal_bam_bai_file': [normal_bam_bai_file],
        'normal_small_variants_vcf_file': [normal_small_variants_vcf_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['savana']['contigs_txt_file'] = savana_contigs_txt_file
    params['savana']['custom_params_file'] = savana_custom_params_file
    params['severus']['vntr_bed_file'] = severus_vntr_bed_file
    params['svisionpro']['model_file'] = svisionpro_model_file
    params['svisionpro']['extra_args'] = '--detect_mode somatic --preset hifi --min_supp 1 --min_mapq 0 --min_sv_size 10 --max_sv_size 1000000 --device cpu'
    params['svisionpro']['extract_extra_args'] = '--extract somatic --min_supp 1'

    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='variant_calling_long-read-dna-somatic.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

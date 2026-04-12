import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


# ClairS, Deepsomatic, Delly2, GRIDSS, Mutect2, Lumpy, Strelka2, Svaba
def test_variant_calling_short_read_dna_somatic_github_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/nexus-dna-002-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    tumor_bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    normal_bam_file = get_data_path(name='bam/nexus-dna-002-normal-paired-end-read_fixmate_markeddup_recalibrated.bam')
    normal_bam_bai_file = get_data_path(name='bam/nexus-dna-002-normal-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/variant_calling_short-read-dna-somatic/params_1.yaml')
    delly2_exclude_tsv_file = get_data_path(name='indices/delly2/human.hg38.excl.tsv')
    octopus_regions_txt_file = get_data_path(name='indices/octopus/hg38_regions.txt')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_shortread_somatic_1'
    work_dir = temp_dir + '/work/test_variant_calling_shortread_somatic_1'
    output_dir = temp_dir + '/outputs/test_variant_calling_shortread_somatic_1'
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
        'tumor_bam_file_no_realign': [tumor_bam_file],
        'tumor_bam_bai_file_no_realign': [tumor_bam_bai_file],
        'normal_bam_file_no_realign': [normal_bam_file],
        'normal_bam_bai_file_no_realign': [normal_bam_bai_file],
        'normal_sample_id': ['nexus-dna-002-normal-paired-end-read']
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['deepsomatic']['input_path'] = '/home/runner/work/nexus/nexus/'
    params['deepsomatic']['output_path'] = '/tmp/'
    params['delly2']['exclude_tsv_file'] = delly2_exclude_tsv_file
    params['mutect2']['germline_resource_vcf_file'] = ''
    params['mutect2']['panel_of_normals_vcf_file'] = ''
    params['mutect2']['getpileupsummaries_variant_vcf_file'] = ''
    params['octopus']['regions_txt_file'] = octopus_regions_txt_file
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='variant_calling_short-read-dna-somatic.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# ClairS, Deepsomatic, Delly2, GRIDSS, Mutect2, Lumpy, Strelka2, Svaba
def test_variant_calling_short_read_dna_somatic_local_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/nexus-dna-002-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    tumor_bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    normal_bam_file = get_data_path(name='bam/nexus-dna-002-normal-paired-end-read_fixmate_markeddup_recalibrated.bam')
    normal_bam_bai_file = get_data_path(name='bam/nexus-dna-002-normal-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/variant_calling_short-read-dna-somatic/params_1.yaml')
    delly2_exclude_tsv_file = get_data_path(name='indices/delly2/human.hg38.excl.tsv')
    octopus_regions_txt_file = get_data_path(name='indices/octopus/hg38_regions.txt')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_shortread_somatic_1'
    work_dir = temp_dir + '/work/test_variant_calling_shortread_somatic_1'
    output_dir = temp_dir + '/outputs/test_variant_calling_shortread_somatic_1'
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
        'tumor_bam_file_no_realign': [tumor_bam_file],
        'tumor_bam_bai_file_no_realign': [tumor_bam_bai_file],
        'normal_bam_file_no_realign': [normal_bam_file],
        'normal_bam_bai_file_no_realign': [normal_bam_bai_file],
        'normal_sample_id': ['nexus-dna-002-normal-paired-end-read']
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['deepsomatic']['input_path'] = '/Users/leework/Documents/Research/projects/project_nexus/'
    params['deepsomatic']['output_path'] = '/var/folders/'
    params['delly2']['exclude_tsv_file'] = delly2_exclude_tsv_file
    params['mutect2']['germline_resource_vcf_file'] = ''
    params['mutect2']['panel_of_normals_vcf_file'] = ''
    params['mutect2']['getpileupsummaries_variant_vcf_file'] = ''
    params['octopus']['regions_txt_file'] = octopus_regions_txt_file
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='variant_calling_short-read-dna-somatic.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# Dysgu
def test_variant_calling_short_read_dna_somatic_2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/nexus-dna-002-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    tumor_bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    normal_bam_file = get_data_path(name='bam/nexus-dna-003-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    normal_bam_bai_file = get_data_path(name='bam/nexus-dna-003-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/variant_calling_short-read-dna-somatic/params_2.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_shortread_somatic_2'
    work_dir = temp_dir + '/work/test_variant_calling_shortread_somatic_2'
    output_dir = temp_dir + '/outputs/test_variant_calling_shortread_somatic_2'
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
        'tumor_bam_file_no_realign': [tumor_bam_file],
        'tumor_bam_bai_file_no_realign': [tumor_bam_bai_file],
        'normal_bam_file_no_realign': [normal_bam_file],
        'normal_bam_bai_file_no_realign': [normal_bam_bai_file],
        'normal_sample_id': ['nexus-dna-002-normal-paired-end-read']
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file

    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='variant_calling_short-read-dna-somatic.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# Manta
def test_variant_calling_short_read_dna_somatic_3():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/nexus-dna-002-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    tumor_bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    normal_bam_file = get_data_path(name='bam/nexus-dna-002-normal-paired-end-read_fixmate_markeddup_recalibrated.bam')
    normal_bam_bai_file = get_data_path(name='bam/nexus-dna-002-normal-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/variant_calling_short-read-dna-somatic/params_3.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_shortread_somatic_3'
    work_dir = temp_dir + '/work/test_variant_calling_shortread_somatic_3'
    output_dir = temp_dir + '/outputs/test_variant_calling_shortread_somatic_3'
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
        'tumor_bam_file_no_realign': [tumor_bam_file],
        'tumor_bam_bai_file_no_realign': [tumor_bam_bai_file],
        'normal_bam_file_no_realign': [normal_bam_file],
        'normal_bam_bai_file_no_realign': [normal_bam_bai_file],
        'normal_sample_id': ['nexus-dna-002-normal-paired-end-read']
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file

    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='variant_calling_short-read-dna-somatic.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# Seqeunza
def test_variant_calling_short_read_dna_somatic_4():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/nexus-dna-005-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    tumor_bam_bai_file = get_data_path(name='bam/nexus-dna-005-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    normal_bam_file = get_data_path(name='bam/nexus-dna-005-normal-paired-end-read_fixmate_markeddup_recalibrated.bam')
    normal_bam_bai_file = get_data_path(name='bam/nexus-dna-005-normal-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr21-22.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/variant_calling_short-read-dna-somatic/params_4.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_shortread_somatic_4'
    work_dir = temp_dir + '/work/test_variant_calling_shortread_somatic_4'
    output_dir = temp_dir + '/outputs/test_variant_calling_shortread_somatic_4'
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
        'tumor_bam_file_no_realign': [tumor_bam_file],
        'tumor_bam_bai_file_no_realign': [tumor_bam_bai_file],
        'normal_bam_file_no_realign': [normal_bam_file],
        'normal_bam_bai_file_no_realign': [normal_bam_bai_file],
        'normal_sample_id': ['nexus-dna-005-normal-paired-end-read']
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file

    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='variant_calling_short-read-dna-somatic.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

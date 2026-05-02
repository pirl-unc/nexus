import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


# DeepVariant + WhatsHap (GitHub-CI variant: deepvariant docker bind-mount paths
# point at /home/runner/work/...).
def test_haplotagging_short_read_dna_github_deepvariant():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_short-read-dna/params_deepvariant.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_shortread_dna_github_deepvariant'
    work_dir = temp_dir + '/work/test_haplotagging_shortread_dna_github_deepvariant'
    output_dir = temp_dir + '/outputs/test_haplotagging_shortread_dna_github_deepvariant'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['deepvariant']['input_path'] = '/home/runner/work/nexus/nexus/'
    params['deepvariant']['output_path'] = '/tmp/'
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='haplotagging_short-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# DeepVariant + WhatsHap (local-laptop variant).
def test_haplotagging_short_read_dna_local_deepvariant():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_short-read-dna/params_deepvariant.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_shortread_dna_local_deepvariant'
    work_dir = temp_dir + '/work/test_haplotagging_shortread_dna_local_deepvariant'
    output_dir = temp_dir + '/outputs/test_haplotagging_shortread_dna_local_deepvariant'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['deepvariant']['input_path'] = '/Users/ajslee/Documents/Research/projects/project_nexus/'
    params['deepvariant']['output_path'] = '/var/folders/'
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='haplotagging_short-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# HaplotypeCaller (chr17 only) + WhatsHap + HapCUT2-WhatsHap.
def test_haplotagging_short_read_dna_haplotypecaller():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_short-read-dna/params_haplotypecaller.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_shortread_dna_haplotypecaller'
    work_dir = temp_dir + '/work/test_haplotagging_shortread_dna_haplotypecaller'
    output_dir = temp_dir + '/outputs/test_haplotagging_shortread_dna_haplotypecaller'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
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
    run_workflow(workflow='haplotagging_short-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# Clair3 (Illumina model) + WhatsHap + HapCUT2-WhatsHap.
def test_haplotagging_short_read_dna_clair3():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_short-read-dna/params_clair3.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_shortread_dna_clair3'
    work_dir = temp_dir + '/work/test_haplotagging_shortread_dna_clair3'
    output_dir = temp_dir + '/outputs/test_haplotagging_shortread_dna_clair3'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
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
    run_workflow(workflow='haplotagging_short-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# Strelka2-germline + WhatsHap + HapCUT2-WhatsHap.
def test_haplotagging_short_read_dna_strelka2_germline():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_short-read-dna/params_strelka2-germline.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_shortread_dna_strelka2_germline'
    work_dir = temp_dir + '/work/test_haplotagging_shortread_dna_strelka2_germline'
    output_dir = temp_dir + '/outputs/test_haplotagging_shortread_dna_strelka2_germline'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
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
    run_workflow(workflow='haplotagging_short-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# All callers (DV + HC + Clair3 + Strelka2-germline) × all phasers
# (WhatsHap + HapCUT2-WhatsHap).
def test_haplotagging_short_read_dna_all_local():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_short-read-dna/params_all.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_shortread_dna_all'
    work_dir = temp_dir + '/work/test_haplotagging_shortread_dna_all'
    output_dir = temp_dir + '/outputs/test_haplotagging_shortread_dna_all'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['deepvariant']['input_path'] = '/Users/ajslee/Documents/Research/projects/project_nexus/'
    params['deepvariant']['output_path'] = '/var/folders/'
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='haplotagging_short-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# All callers (DV + HC + Clair3 + Strelka2-germline) × all phasers
# (WhatsHap + HapCUT2-WhatsHap), GitHub-CI variant: deepvariant docker
# bind-mount paths point at /home/runner/work/...
def test_haplotagging_short_read_dna_all_github():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_short-read-dna/params_all.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_shortread_dna_all_github'
    work_dir = temp_dir + '/work/test_haplotagging_shortread_dna_all_github'
    output_dir = temp_dir + '/outputs/test_haplotagging_shortread_dna_all_github'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['deepvariant']['input_path'] = '/home/runner/work/nexus/nexus/'
    params['deepvariant']['output_path'] = '/tmp/'
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='haplotagging_short-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


# Clair3, Longshot, NanoCaller, NanoVar (SNV callers — nexus-dna-001-tumor)
def test_variant_calling_long_read_dna_germline_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/variant_calling_long-read-dna-germline/params_1.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_longread_germline_1'
    work_dir = temp_dir + '/work/test_variant_calling_longread_germline_1'
    output_dir = temp_dir + '/outputs/test_variant_calling_longread_germline_1'
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
    run_workflow(workflow='variant_calling_long-read-dna-germline.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# CuteSV, Delly2, Dysgu, PBSV, Sniffles2, SVIM (SV callers — nexus-dna-004-tumor)
def test_variant_calling_long_read_dna_germline_2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-004-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-004-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/variant_calling_long-read-dna-germline/params_2.yaml')
    delly2_exclude_tsv_file = get_data_path(name='indices/delly2/human.hg38.excl.tsv')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_longread_germline_2'
    work_dir = temp_dir + '/work/test_variant_calling_longread_germline_2'
    output_dir = temp_dir + '/outputs/test_variant_calling_longread_germline_2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-004-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
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
    run_workflow(workflow='variant_calling_long-read-dna-germline.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# DeepVariant + HiFiCNV (nexus-dna-001-tumor)
def test_variant_calling_long_read_dna_germline_github_3():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-005-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-005-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr21-22.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/variant_calling_long-read-dna-germline/params_3.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_longread_germline_3'
    work_dir = temp_dir + '/work/test_variant_calling_longread_germline_3'
    output_dir = temp_dir + '/outputs/test_variant_calling_longread_germline_3'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-005-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    # Create a minimal exclude BED for HiFiCNV
    exclude_bed_file = intermediate_dir + '/hificnv_exclude.bed'
    with open(exclude_bed_file, 'w') as f:
        f.write("")

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['deepvariant']['input_path'] = '/home/runner/work/nexus/nexus/'
    params['deepvariant']['output_path'] = '/tmp/'
    params['hificnv']['exclude_bed_file'] = exclude_bed_file

    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='variant_calling_long-read-dna-germline.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


# DeepVariant + HiFiCNV (nexus-dna-001-tumor)
def test_variant_calling_long_read_dna_germline_local_3():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-005-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-005-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr21-22.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/variant_calling_long-read-dna-germline/params_3.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_longread_germline_3'
    work_dir = temp_dir + '/work/test_variant_calling_longread_germline_3'
    output_dir = temp_dir + '/outputs/test_variant_calling_longread_germline_3'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-005-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    # Create a minimal exclude BED for HiFiCNV
    exclude_bed_file = intermediate_dir + '/hificnv_exclude.bed'
    with open(exclude_bed_file, 'w') as f:
        f.write("")

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params['reference_genome_fasta_file'] = reference_genome_fasta_file
    params['deepvariant']['input_path'] = '/Users/leework/Documents/Research/projects/project_nexus/'
    params['deepvariant']['output_path'] = '/var/folders/'
    params['hificnv']['exclude_bed_file'] = exclude_bed_file

    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='variant_calling_long-read-dna-germline.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

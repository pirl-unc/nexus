import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_haplotype_phasing_long_read_dna_github_deepvariant():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotype_phasing_long-read-dna/params_deepvariant.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotype_phasing_longread_dna_github_deepvariant'
    work_dir = temp_dir + '/work/test_haplotype_phasing_longread_dna_github_deepvariant'
    output_dir = temp_dir + '/outputs/test_haplotype_phasing_longread_dna_github_deepvariant'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # sample_id must match the BAM's @RG SM tag because hiphase validates
    # --sample-name against the VCF sample column (which DeepVariant pulls
    # from the BAM's read group). The BAM's SM tag is 'nexus-dna-002-tumor-long-read'.
    pd.DataFrame({
        'sample_id': ['nexus-dna-002-tumor'],
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
    run_workflow(workflow='haplotype_phasing_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_haplotype_phasing_long_read_dna_local_deepvariant():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotype_phasing_long-read-dna/params_deepvariant.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotype_phasing_longread_dna_local_deepvariant'
    work_dir = temp_dir + '/work/test_haplotype_phasing_longread_dna_local_deepvariant'
    output_dir = temp_dir + '/outputs/test_haplotype_phasing_longread_dna_local_deepvariant'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-002-tumor-long-read'],
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
    run_workflow(workflow='haplotype_phasing_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_haplotype_phasing_long_read_dna_longshot():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotype_phasing_long-read-dna/params_longshot.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotype_phasing_longread_dna_longshot'
    work_dir = temp_dir + '/work/test_haplotype_phasing_longread_dna_longshot'
    output_dir = temp_dir + '/outputs/test_haplotype_phasing_longread_dna_longshot'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # sample_id must match the BAM's @RG SM tag for hiphase --sample-name validation.
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor-long-read'],
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
    run_workflow(workflow='haplotype_phasing_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_haplotype_phasing_long_read_dna_all_github():
    """Run all small-variants callers (DeepVariant + Longshot) AND all phaser
    methods (pbsv + hiphase + whatshap + hapcut2-whatshap + standalone longshot).

    Caller-dependent phasers (hiphase, whatshap, hapcut2-whatshap) run twice —
    once per caller's VCF — with outputs namespaced under
    with_deepvariant/ and with_longshot/ subdirs of output_dir.

    Uses dna-001 because Longshot needs enough variants on chr17 to produce
    a non-empty VCF for HiPhase to consume.
    """
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotype_phasing_long-read-dna/params_all.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotype_phasing_longread_dna_all'
    work_dir = temp_dir + '/work/test_haplotype_phasing_longread_dna_all'
    output_dir = temp_dir + '/outputs/test_haplotype_phasing_longread_dna_all'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # sample_id must match the BAM's @RG SM tag for hiphase --sample-name validation.
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor-long-read'],
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
    run_workflow(workflow='haplotype_phasing_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_haplotype_phasing_long_read_dna_all_local():
    """Run all small-variants callers (DeepVariant + Longshot) AND all phaser
    methods (pbsv + hiphase + whatshap + hapcut2-whatshap + standalone longshot).

    Caller-dependent phasers (hiphase, whatshap, hapcut2-whatshap) run twice —
    once per caller's VCF — with outputs namespaced under
    with_deepvariant/ and with_longshot/ subdirs of output_dir.

    Uses dna-001 because Longshot needs enough variants on chr17 to produce
    a non-empty VCF for HiPhase to consume.
    """
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotype_phasing_long-read-dna/params_all.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotype_phasing_longread_dna_all'
    work_dir = temp_dir + '/work/test_haplotype_phasing_longread_dna_all'
    output_dir = temp_dir + '/outputs/test_haplotype_phasing_longread_dna_all'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # sample_id must match the BAM's @RG SM tag for hiphase --sample-name validation.
    pd.DataFrame({
        'sample_id': ['nexus-dna-002-tumor-long-read'],
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
    run_workflow(workflow='haplotype_phasing_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_haplotype_phasing_long_read_dna_hapcut2_whatshap_local():
    """Focused test for the HapCUT2-WhatsHap phaser path.

    Runs only the hapcut2-whatshap method (no hiphase, no whatshap, no pbsv,
    no standalone longshot). Uses DeepVariant as the small-variants caller.
    Verifies that:
      1. extractHAIRS + HAPCUT2 produce a phased VCF.
      2. WhatsHap haplotag produces a haplotagged BAM using HapCUT2's VCF.
    """
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotype_phasing_long-read-dna/params_hapcut2-whatshap.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotype_phasing_longread_dna_hapcut2_whatshap'
    work_dir = temp_dir + '/work/test_haplotype_phasing_longread_dna_hapcut2_whatshap'
    output_dir = temp_dir + '/outputs/test_haplotype_phasing_longread_dna_hapcut2_whatshap'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor-long-read'],
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
    run_workflow(workflow='haplotype_phasing_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_haplotype_phasing_long_read_dna_longshot_phaser_only():
    """Focused test for the standalone Longshot phaser path.

    Verifies that requesting longshot via the `methods` list (with DeepVariant
    as the small-variants caller) makes Longshot run and publish its phased
    VCF + BAM standalone — without chaining its VCF into HiPhase / WhatsHap /
    HapCUT2-WhatsHap.

    Pipeline that runs:
      - DeepVariant (caller; its VCF is unused since no caller-dependent phaser
        is selected by methods)
      - Longshot   (standalone phaser; publishes phased VCF + BAM)
    """
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotype_phasing_long-read-dna/params_longshot.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotype_phasing_longread_dna_longshot_phaser_only'
    work_dir = temp_dir + '/work/test_haplotype_phasing_longread_dna_longshot_phaser_only'
    output_dir = temp_dir + '/outputs/test_haplotype_phasing_longread_dna_longshot_phaser_only'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor-long-read'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    # Use the Longshot-only fixture: Longshot is the small-variants caller and
    # publishes its phased VCF + BAM directly, with no additional phaser methods.
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
    run_workflow(workflow='haplotype_phasing_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

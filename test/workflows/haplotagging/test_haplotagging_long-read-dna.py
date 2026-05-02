import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_haplotagging_long_read_dna_github_deepvariant():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_long-read-dna/params_deepvariant.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_longread_dna_github_deepvariant'
    work_dir = temp_dir + '/work/test_haplotagging_longread_dna_github_deepvariant'
    output_dir = temp_dir + '/outputs/test_haplotagging_longread_dna_github_deepvariant'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
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
    run_workflow(workflow='haplotagging_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_haplotagging_long_read_dna_local_deepvariant():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_long-read-dna/params_deepvariant.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_longread_dna_local_deepvariant'
    work_dir = temp_dir + '/work/test_haplotagging_longread_dna_local_deepvariant'
    output_dir = temp_dir + '/outputs/test_haplotagging_longread_dna_local_deepvariant'
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
    run_workflow(workflow='haplotagging_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_haplotagging_long_read_dna_longshot():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_long-read-dna/params_longshot.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_longread_dna_longshot'
    work_dir = temp_dir + '/work/test_haplotagging_longread_dna_longshot'
    output_dir = temp_dir + '/outputs/test_haplotagging_longread_dna_longshot'
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
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='haplotagging_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_haplotagging_long_read_dna_all_github():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_long-read-dna/params_all.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_longread_dna_all'
    work_dir = temp_dir + '/work/test_haplotagging_longread_dna_all'
    output_dir = temp_dir + '/outputs/test_haplotagging_longread_dna_all'
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
    run_workflow(workflow='haplotagging_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_haplotagging_long_read_dna_all_local():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    margin_phase_params_json_file = get_data_path(
        name='indices/margin/phase/allParams.phase_vcf.pb-hifi.json'
    )
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_long-read-dna/params_all.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_longread_dna_all'
    work_dir = temp_dir + '/work/test_haplotagging_longread_dna_all'
    output_dir = temp_dir + '/outputs/test_haplotagging_longread_dna_all'
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
    params.setdefault('margin', {})['phase_params_json_file'] = margin_phase_params_json_file
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='haplotagging_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_haplotagging_long_read_dna_all_tsv_local():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    margin_phase_params_json_file = get_data_path(
        name='indices/margin/phase/allParams.phase_vcf.pb-hifi.json'
    )
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_long-read-dna/params_all_tsv.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_longread_dna_all_tsv'
    work_dir = temp_dir + '/work/test_haplotagging_longread_dna_all_tsv'
    output_dir = temp_dir + '/outputs/test_haplotagging_longread_dna_all_tsv'
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
    params.setdefault('margin', {})['phase_params_json_file'] = margin_phase_params_json_file
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='haplotagging_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

    sample_out = os.path.join(output_dir, 'nexus-dna-002-tumor-long-read')

    leaked_bams = []
    for root, _, fnames in os.walk(sample_out):
        for fn in fnames:
            if fn.endswith('_haplotagged.bam') or fn.endswith('_phased.bam') \
               or fn.endswith('_haplotagged.bam.bai') or fn.endswith('_phased.bam.bai'):
                leaked_bams.append(os.path.join(root, fn))
    assert not leaked_bams, \
        f"haplotag_output='tsv' but found published haplotagged BAMs: {leaked_bams}"

    expected_tsvs = [
        os.path.join(sample_out, 'deepvariant_whatshap'),
        os.path.join(sample_out, 'longshot_whatshap'),
        os.path.join(sample_out, 'clair3_whatshap'),
        os.path.join(sample_out, 'deepvariant_pbsv_hiphase'),
        os.path.join(sample_out, 'longshot_pbsv_hiphase'),
        os.path.join(sample_out, 'longshot_hapcut2-whatshap'),
        os.path.join(sample_out, 'clair3_hapcut2-whatshap'),
        os.path.join(sample_out, 'deepvariant_longphase'),
        os.path.join(sample_out, 'longshot_longphase'),
        os.path.join(sample_out, 'clair3_longphase'),
    ]
    for d in expected_tsvs:
        assert os.path.isdir(d), f"missing phaser output dir: {d}"
        tsvs = [fn for fn in os.listdir(d) if fn.endswith('.tsv') or fn.endswith('.tsv.gz')]
        assert tsvs, f"no haplotag TSV found in {d} (contents: {os.listdir(d)})"


def test_haplotagging_long_read_dna_hapcut2_whatshap_local():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_long-read-dna/params_hapcut2-whatshap.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_longread_dna_hapcut2_whatshap'
    work_dir = temp_dir + '/work/test_haplotagging_longread_dna_hapcut2_whatshap'
    output_dir = temp_dir + '/outputs/test_haplotagging_longread_dna_hapcut2_whatshap'
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
    run_workflow(workflow='haplotagging_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_haplotagging_long_read_dna_longphase_local():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_long-read-dna/params_longphase.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_longread_dna_longphase'
    work_dir = temp_dir + '/work/test_haplotagging_longread_dna_longphase'
    output_dir = temp_dir + '/outputs/test_haplotagging_longread_dna_longphase'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    sample_id = 'nexus-dna-001-tumor-long-read'
    pd.DataFrame({
        'sample_id': [sample_id],
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
    run_workflow(workflow='haplotagging_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

    longphase_dir = os.path.join(output_dir, sample_id, 'deepvariant_longphase')
    phased_vcf = os.path.join(longphase_dir, f'{sample_id}_longphase_phased.vcf.gz')
    phased_vcf_tbi = phased_vcf + '.tbi'
    bam_basename = os.path.basename(bam_file).replace('.bam', '')
    haplotagged_bam = os.path.join(longphase_dir, f'{bam_basename}_haplotagged.bam')
    haplotagged_bai = haplotagged_bam + '.bai'
    assert os.path.isdir(longphase_dir), f"Expected LongPhase output dir at {longphase_dir}"
    assert os.path.exists(phased_vcf), f"Expected phased VCF at {phased_vcf}"
    assert os.path.exists(phased_vcf_tbi), f"Expected phased VCF index at {phased_vcf_tbi}"
    assert os.path.exists(haplotagged_bam), f"Expected haplotagged BAM at {haplotagged_bam}"
    assert os.path.exists(haplotagged_bai), f"Expected haplotagged BAM index at {haplotagged_bai}"


def test_haplotagging_long_read_dna_longshot_phaser_only():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/haplotagging_long-read-dna/params_longshot.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_longread_dna_longshot_phaser_only'
    work_dir = temp_dir + '/work/test_haplotagging_longread_dna_longshot_phaser_only'
    output_dir = temp_dir + '/outputs/test_haplotagging_longread_dna_longshot_phaser_only'
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
    run_workflow(workflow='haplotagging_long-read-dna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

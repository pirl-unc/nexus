import gzip
import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


# GeLuster, isONclust3 (sample)
def test_read_clustering_longread_rna():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fastq_file = get_data_path(name='inputs/workflows/assembly_long-read-rna/reads.fastq.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/read_clustering_long-read-rna/params_1.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_read_clustering_longread_rna'
    work_dir = temp_dir + '/work/test_read_clustering_longread_rna'
    output_dir = temp_dir + '/outputs/test_read_clustering_longread_rna'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample'],
        'fastq_file': [fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file'] = f"{intermediate_dir}/samples.tsv"
    params['output_dir'] = output_dir
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file,
    ]
    run_workflow(workflow='read_clustering_long-read-rna.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

    # Every method must publish a normalized per-read cluster TSV that assigns
    # each input read to exactly one cluster. This is the contract the
    # assembly_*_clustered subworkflows rely on, and a workflow that merely
    # completes does not prove it. GeLuster in particular reports singleton
    # reads in a second file, so a naive parse silently drops them.
    read_names = set()
    with gzip.open(fastq_file, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 0 and len(line) > 1:
                read_names.add(line[1:].split()[0])
    for method in ['geluster', 'isonclust3']:
        tsv_file = f"{output_dir}/sample/{method}/sample_{method}_clusters.tsv"
        assert os.path.exists(tsv_file), f"{method} did not publish {tsv_file}"
        df_clusters = pd.read_csv(tsv_file, sep='\t')
        assert list(df_clusters.columns) == ['cluster_id', 'read_name']
        clustered = df_clusters['read_name'].astype(str)
        assert not clustered.duplicated().any(), f"{method} assigned a read twice"
        unknown = set(clustered) - read_names
        assert not unknown, f"{method} emitted unknown read names: {sorted(unknown)[:5]}"
        missing = read_names - set(clustered)
        assert not missing, f"{method} dropped {len(missing)} reads, e.g. {sorted(missing)[:5]}"

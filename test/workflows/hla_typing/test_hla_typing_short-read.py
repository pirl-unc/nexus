import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


# arcasHLA, HLAProfiler, seq2HLA (nexus-rna-104-normal)
def test_hla_typing_shortread():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-rna-104-normal-paired-end-read_star_outputs/nexus-rna-104-normal-paired-end-read_star_Aligned.sortedByCoord.out.bam')
    bam_bai_file = get_data_path(name='bam/nexus-rna-104-normal-paired-end-read_star_outputs/nexus-rna-104-normal-paired-end-read_star_Aligned.sortedByCoord.out.bam.bai')
    fastq_file_1 = get_data_path(name='fastq/nexus-rna-104-normal_paired-end_read_r1.fastq.gz')
    fastq_file_2 = get_data_path(name='fastq/nexus-rna-104-normal_paired-end_read_r2.fastq.gz')
    params_yaml_file = get_data_path(name='inputs/workflows/hla_typing_short-read/params_1.yaml')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_hla_typing_shortread_local_1'
    work_dir = temp_dir + '/work/test_hla_typing_shortread_local_1'
    output_dir = temp_dir + '/outputs/test_hla_typing_shortread_local_1'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-104-normal'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file],
        'fastq_file_1': [fastq_file_1],
        'fastq_file_2': [fastq_file_2]
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
    run_workflow(workflow='hla_typing_short-read.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

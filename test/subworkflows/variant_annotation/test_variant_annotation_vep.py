import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_annotation_vep():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    vcf_file = get_data_path(name='vcf/nexus-dna-001-tumor_gatk4-mutect2.vcf')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_annotation_vep'
    work_dir = temp_dir + '/work/test_variant_annotation_vep'
    output_dir = temp_dir + '/outputs/test_variant_annotation_vep'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor'],
        'vcf_file': [vcf_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--vep_dir', '/Users/leework/.vep',
        '--params_vep', '"--species homo_sapiens --database --offline --cache --assembly GRCh38"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_annotation_vep.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

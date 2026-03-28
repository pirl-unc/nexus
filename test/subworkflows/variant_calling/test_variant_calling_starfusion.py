import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_calling_starfusion():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    fastq_file_1 = get_data_path(name='fastq/nexus-rna-004-tumor_paired-end_read_r1.fastq.gz')
    fastq_file_2 = get_data_path(name='fastq/nexus-rna-004-tumor_paired-end_read_r2.fastq.gz')
    genome_lib_dir = '/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/starfusion/ctat_genome_lib_starfv1.10/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/'
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_starfusion'
    work_dir = temp_dir + '/work/test_variant_calling_starfusion'
    output_dir = temp_dir + '/outputs/test_variant_calling_starfusion'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-004-tumor'],
        'fastq_file_1': [fastq_file_1],
        'fastq_file_2': [fastq_file_2],
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--genome_lib_dir', genome_lib_dir,
        '--params_starfusion', '"--full_Monty --no_filter --STAR_twopass"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_starfusion.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

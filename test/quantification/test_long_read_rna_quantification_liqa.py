import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path, get_alias_path


def test_long_read_rna_quantification_liqa():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/sample200tumor_minimap2_mdtagged_sorted.bam')
    tumor_bam_bai_file = get_data_path(name='bam/sample200tumor_minimap2_mdtagged_sorted.bam.bai')
    normal_bam_file = get_data_path(name='bam/sample200normal_minimap2_mdtagged_sorted.bam')
    normal_bam_bai_file = get_data_path(name='bam/sample200normal_minimap2_mdtagged_sorted.bam.bai')
    liqa_refgene_file = get_data_path(name='indices/liqa/liqa_gencode_v41_tp53.refgene')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_rna_quantification_liqa'
    work_dir = temp_dir + '/work/test_long_read_rna_quantification_liqa'
    output_dir = temp_dir + '/outputs/test_long_read_rna_quantification_liqa'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample200tumor', 'sample200normal'],
        'bam_file': [tumor_bam_file, normal_bam_file],
        'bam_bai_file': [tumor_bam_bai_file, normal_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--liqa_refgene_file', liqa_refgene_file,
        '--params_liqa_quantify', '"-max_distance 20 -f_weight 1"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='long_read_rna_quantification_liqa.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

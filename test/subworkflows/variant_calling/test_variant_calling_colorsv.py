import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_calling_colorsv():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    gfa_file = get_data_path(name='gfa/colorsv_demo.asm.bp.r_utg.gfa')
    bed_file = get_data_path(name='bed/colorsv_chm13v2.cen-mask.bed')
    reference_genome_fasta_file = "/Users/leework/Documents/Research/projects/seqdata/references/chm13v2.0.fa"
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_colorsv'
    work_dir = temp_dir + '/work/test_variant_calling_colorsv'
    output_dir = temp_dir + '/outputs/test_variant_calling_colorsv'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample'],
        'gfa_file': [gfa_file],
        'tumor_ids': ['m84039_230312_025934_s1 m84039_230328_000836_s3']
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--filter_bed_file', bed_file,
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_colorsv.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

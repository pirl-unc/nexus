import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_transcriptome_assembly_rnabloom2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_rna_fastq_file = get_data_path(name='fastq/sample200tumor_long_read_rna.fastq.gz')
    normal_rna_fastq_file = get_data_path(name='fastq/sample200normal_long_read_rna.fastq.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_transcriptome_assembly_rnabloom2'
    work_dir = temp_dir + '/work/test_transcriptome_assembly_rnabloom2'
    output_dir = temp_dir + '/outputs/test_transcriptome_assembly_rnabloom2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample200tumor', 'sample200normal'],
        'fastq_file': [tumor_rna_fastq_file, normal_rna_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--params_rnabloom2', '"--qual 20 --qual-avg 20 --mincov 3 -ntcard -savebf"',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='transcriptome_assembly_rnabloom2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

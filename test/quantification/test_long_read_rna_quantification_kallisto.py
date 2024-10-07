import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path, get_alias_path


def test_long_read_rna_quantification_kallisto():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_fastq_file = get_data_path(name='fastq/sample200tumor_long_read_rna.fastq.gz')
    normal_fastq_file = get_data_path(name='fastq/sample200normal_long_read_rna.fastq.gz')
    kallisto_index_file = get_data_path(name='indices/kallisto/kallisto_gencode_v41_tp53_index_k63.idx')
    gtf_file = get_data_path(name='gtf/gencode_v41_tp53_annotation.gtf')
    t2g_file = get_data_path(name='indices/kallisto/gencode_v41.t2g')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_long_read_rna_quantification_kallisto'
    work_dir = temp_dir + '/work/test_long_read_rna_quantification_kallisto'
    output_dir = temp_dir + '/outputs/test_long_read_rna_quantification_kallisto'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample200tumor','sample200normal'],
        'fastq_file': [tumor_fastq_file, normal_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--kallisto_index_file', kallisto_index_file,
        '--gtf_file', gtf_file,
        '--t2g_file', t2g_file,
        '--params_kallisto_bus', '"-x bulk --threshold 0.8"',
        '--params_bustools_sort', '""',
        '--params_bustools_count', '"--cm -m"',
        '--params_kallisto_quanttcc', '"-P PacBio --matrix-to-files"',
        '--kallisto_index_file', kallisto_index_file,
        '--output_dir', output_dir
    ]
    run_workflow(workflow='long_read_rna_quantification_kallisto.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

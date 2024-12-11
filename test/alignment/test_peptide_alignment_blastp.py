import pandas as pd
import os
from nexuslib.main import run_workflow
from ..data import get_data_path


def test_peptide_alignment_blastp():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    peptides_fasta_file = get_data_path(name='fasta/hg38_tp53_sampled_peptides.fa')
    reference_proteome_fasta_file = get_data_path(name='fasta/gencode_v41_tp53_translations.fa')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_peptide_alignment_blastp'
    work_dir = temp_dir + '/work/test_peptide_alignment_blastp'
    output_dir = temp_dir + '/outputs/test_peptide_alignment_blastp'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['sample300tumor'],
        'fasta_file': [peptides_fasta_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_proteome_fasta_file', reference_proteome_fasta_file,
        '--output_dir', output_dir
    ]
    run_workflow(workflow='peptide_alignment_blastp.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


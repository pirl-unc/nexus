import pandas as pd
import os
import yaml
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_peptide_prediction_exacto():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')

    # ---- Long-read FASTQ.GZ inputs (one tumor DNA, one matched-normal DNA,
    #      one tumor RNA, all from the same chr17 simulated sample) ----
    tumor_dna_fastq_file  = get_data_path(name='fastq/nexus-dna-001-tumor_long_read.fastq.gz')
    normal_dna_fastq_file = get_data_path(name='fastq/nexus-dna-001-normal_long_read.fastq.gz')
    tumor_rna_fastq_file  = get_data_path(name='fastq/nexus-rna-002-tumor_long_read.fastq.gz')

    # ---- Reference files (chr17 subset shared with the mopepgen test) ----
    reference_genome_fasta_file    = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    reference_gene_annotation_file = get_data_path(name='gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz')
    reference_proteome_fasta_file  = get_data_path(name='fasta/gencode.v45.pc_translations.fa')

    # ---- Params template ----
    params_yaml_file = get_data_path(name='inputs/subworkflows/peptide_prediction_exacto/params.yaml')

    # ---- Test scratch dirs ----
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_peptide_prediction_exacto'
    work_dir         = temp_dir + '/work/test_peptide_prediction_exacto'
    output_dir       = temp_dir + '/outputs/test_peptide_prediction_exacto'
    for d in (intermediate_dir, work_dir, output_dir):
        if not os.path.exists(d):
            os.makedirs(d)

    # ---- Build samples TSV ----
    pd.DataFrame({
        'sample_id':             ['nexus-001'],
        'tumor_dna_fastq_file':  [tumor_dna_fastq_file],
        'normal_dna_fastq_file': [normal_dna_fastq_file],
        'tumor_rna_fastq_file':  [tumor_rna_fastq_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)

    # ---- Materialize params.yaml with test-specific paths ----
    with open(params_yaml_file, 'r') as f:
        params = yaml.safe_load(f)
    params['samples_tsv_file']               = f"{intermediate_dir}/samples.tsv"
    params['output_dir']                     = output_dir
    params['reference_genome_fasta_file']    = reference_genome_fasta_file
    params['reference_gene_annotation_file'] = reference_gene_annotation_file
    params['reference_proteome_fasta_file']  = reference_proteome_fasta_file
    params_file = intermediate_dir + '/params.yaml'
    with open(params_file, 'w') as f:
        yaml.dump(params, f, default_flow_style=False, default_style='"')

    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '-params-file', params_file
    ]
    run_workflow(workflow='peptide_prediction_exacto.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

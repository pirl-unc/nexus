import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_calling_sequenza():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    tumor_bam_file = get_data_path(name='bam/nexus-dna-005-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    tumor_bam_bai_file = get_data_path(name='bam/nexus-dna-005-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    normal_bam_file = get_data_path(name='bam/nexus-dna-005-normal-paired-end-read_fixmate_markeddup_recalibrated.bam')
    normal_bam_bai_file = get_data_path(name='bam/nexus-dna-005-normal-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr21-22.fa.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_sequenza'
    work_dir = temp_dir + '/work/test_variant_calling_sequenza'
    output_dir = temp_dir + '/outputs/test_variant_calling_sequenza'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-005-tumor'],
        'tumor_bam_file': [tumor_bam_file],
        'tumor_bam_bai_file': [tumor_bam_bai_file],
        'normal_bam_file': [normal_bam_file],
        'normal_bam_bai_file': [normal_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--assembly', '"hg38"',
        '--chromosomes', '"chr21 chr22"',
        '--sequenzautils_gcwiggle', '"-w 50"',
        '--sequenzautils_bam2seqz', '"-N 1 --qformat sanger --hom 0.8 --het 0.4 --qlimit 1"',
        '--sequenzautils_seqzbinning', '"--window 50"',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_sequenza.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

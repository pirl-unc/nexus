import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_variant_calling_delly2_sr_germline():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-002-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-002-tumor-paired-end-read_fixmate_markeddup_recalibrated.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    exclude_tsv_file = get_data_path(name='indices/delly2/human.hg38.excl.tsv')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_variant_calling_delly2_sr_germline'
    work_dir = temp_dir + '/work/test_variant_calling_delly2_sr_germline'
    output_dir = temp_dir + '/outputs/test_variant_calling_delly2_sr_germline'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-002-tumor'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--exclude_tsv_file', exclude_tsv_file,
        '--params_delly2call', '""',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='variant_calling_delly2-sr-germline.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

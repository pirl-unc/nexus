import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_haplotagging_longshot():
    """End-to-end test for the standalone Longshot haplotagging subworkflow.

    Longshot calls AND phases small variants in a single step, producing both
    a phased VCF and a haplotagged BAM. Uses dna-001 because dna-002 has too
    few variants on chr17 for Longshot to call any at default --min_cov.
    """
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    bam_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam')
    bam_bai_file = get_data_path(name='bam/nexus-dna-001-tumor-long-read_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_haplotagging_longshot'
    work_dir = temp_dir + '/work/test_haplotagging_longshot'
    output_dir = temp_dir + '/outputs/test_haplotagging_longshot'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor-long-read'],
        'bam_file': [bam_file],
        'bam_bai_file': [bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='haplotagging_longshot.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

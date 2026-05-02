import pandas as pd
import os
from nexuslib.main import run_workflow
from ...data import get_data_path


def test_alignment_minimap2_1():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_dna_fastq_file_1 = get_data_path(name='fastq/nexus-dna-001-tumor_long_read.fastq.gz')
    long_read_dna_fastq_file_2 = get_data_path(name='fastq/nexus-dna-002-tumor_long_read.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_alignment_minimap2_1'
    work_dir = temp_dir + '/work/test_alignment_minimap2_1'
    output_dir = temp_dir + '/outputs/test_alignment_minimap2_1'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-dna-001-tumor', 'nexus-dna-002-tumor'],
        'fastq_file': [long_read_dna_fastq_file_1, long_read_dna_fastq_file_2]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--params_minimap2', '"-ax map-hifi --ds --cs --eqx -Y -L --secondary=no"',
        '--platform_tag', 'pacbio',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='alignment_minimap2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_alignment_minimap2_2():
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_rna_fastq_file_1 = get_data_path(name='fastq/nexus-rna-001-tumor_long_read.fastq.gz')
    long_read_rna_fastq_file_2 = get_data_path(name='fastq/nexus-rna-002-tumor_long_read.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_alignment_minimap2_2'
    work_dir = temp_dir + '/work/test_alignment_minimap2_2'
    output_dir = temp_dir + '/outputs/test_alignment_minimap2_2'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['nexus-rna-001-tumor', 'nexus-rna-002-tumor'],
        'fastq_file': [long_read_rna_fastq_file_1, long_read_rna_fastq_file_2]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--params_minimap2', '"-ax splice:hq -uf --cs --eqx -Y -L --secondary=no"',
        '--platform_tag', 'pacbio',
        '--output_dir', output_dir
    ]
    run_workflow(workflow='alignment_minimap2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)


def test_alignment_minimap2_multi_fastq_per_sample():
    """Two TSV rows share sample_id='nexus-dna-001-tumor' but each row points
    to a different fastq file. Both fastqs should be aligned together into a
    single merged BAM via groupTuple — exercising the multi-row-per-sample-id
    code path. A second sample (single row) verifies the original
    one-fastq-per-sample case still works alongside.
    """
    nextflow_config_file = get_data_path(name='nextflow/nextflow_test_docker.config')
    long_read_dna_fastq_file_1 = get_data_path(name='fastq/nexus-dna-001-tumor_long_read.fastq.gz')
    long_read_dna_fastq_file_2 = get_data_path(name='fastq/nexus-dna-002-tumor_long_read.fastq.gz')
    reference_genome_fasta_file = get_data_path(name='fasta/GRCh38.p14.genome.chr17.fa.gz')
    temp_dir = os.getcwd() + '/tmp'
    intermediate_dir = temp_dir + '/intermediate/test_alignment_minimap2_multi_fastq'
    work_dir = temp_dir + '/work/test_alignment_minimap2_multi_fastq'
    output_dir = temp_dir + '/outputs/test_alignment_minimap2_multi_fastq'
    for d in (intermediate_dir, work_dir, output_dir):
        if not os.path.exists(d):
            os.makedirs(d)
    # Two rows for 'nexus-dna-merged' (different fastq files, same sample_id),
    # plus a single row for 'nexus-dna-002-tumor' as a control.
    pd.DataFrame({
        'sample_id': [
            'nexus-dna-merged',
            'nexus-dna-merged',
            'nexus-dna-002-tumor',
        ],
        'fastq_file': [
            long_read_dna_fastq_file_1,
            long_read_dna_fastq_file_2,
            long_read_dna_fastq_file_2,
        ]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '-w', work_dir,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--params_minimap2', '"-ax map-hifi --cs --eqx -Y -L --secondary=no"',
        '--platform_tag', 'pacbio',
        '--output_dir', output_dir,
    ]
    run_workflow(workflow='alignment_minimap2.nf',
                 nextflow='nextflow',
                 workflow_args=workflow_args)

    # Exactly one merged BAM should appear for the multi-row sample_id, and
    # one BAM for the single-row control. (If groupTuple failed, we'd get
    # two BAMs for the merged sample, or output collisions.)
    merged_bam = os.path.join(output_dir, 'nexus-dna-merged_minimap2_sorted.bam')
    merged_bai = merged_bam + '.bai'
    control_bam = os.path.join(output_dir, 'nexus-dna-002-tumor_minimap2_sorted.bam')
    control_bai = control_bam + '.bai'
    assert os.path.exists(merged_bam), \
        f"Expected merged BAM at {merged_bam}"
    assert os.path.exists(merged_bai), \
        f"Expected merged BAM index at {merged_bai}"
    assert os.path.exists(control_bam), \
        f"Expected control BAM at {control_bam}"
    assert os.path.exists(control_bai), \
        f"Expected control BAM index at {control_bai}"

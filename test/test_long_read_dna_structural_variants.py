import pandas as pd
import os
from .data import get_data_path
from nexuslib.main import run_workflow


def test_long_read_dna_structural_variants():
    nextflow_config_file = get_data_path(name='nextflow_test.config')
    long_read_tumor_dna_bam_file = get_data_path(name='hg38_tp53_tumor_dna_minimap2_mdtagged_sorted.bam')
    long_read_tumor_dna_bam_bai_file = get_data_path(name='hg38_tp53_tumor_dna_minimap2_mdtagged_sorted.bam.bai')
    reference_genome_fasta_file = get_data_path(name='hg38_chr17_1-8000000.fa')
    temp_dir = os.getcwd() + '/tmp/'
    intermediate_dir = temp_dir + 'intermediate/test_long_read_dna_structural_variants/'
    work_dir = temp_dir + 'work/test_long_read_dna_structural_variants/'
    output_dir = temp_dir + 'outputs/test_long_read_dna_structural_variants/'
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pd.DataFrame({
        'sample_id': ['tumor'],
        'bam_file': [long_read_tumor_dna_bam_file],
        'bam_bai_file': [long_read_tumor_dna_bam_bai_file]
    }).to_csv(intermediate_dir + "/samples.tsv", sep='\t', index=False)
    workflow_args = [
        '-c', nextflow_config_file,
        '--samples_tsv_file', intermediate_dir + '/samples.tsv',
        '--reference_genome_fasta_file', reference_genome_fasta_file,
        '--tools_list', 'sniffles2,svim,cutesv',    # pbsv is only available on linux
        '--sniffles2', 'sniffles',
        '--sniffles2_params', '"--minsupport 3 --minsvlen 30 --mapq 20 --output-rnames "',
        '--pbsv', 'pbsv',
        '--pbsv_discover_params', '"--ccs --min-gap-comp-id-perc 97 --min-mapq 20 "',
        '--pbsv_call_params', '"--ccs --call-min-reads-per-strand-all-samples 0 --call-min-read-perc-one-sample 10 --call-min-reads-all-samples 3 --call-min-reads-one-sample 3 "',
        '--svim', 'svim',
        '--svim_params', '"--min_mapq 20 --min_sv_size 30 --insertion_sequences --read_names --zmws "',
        '--cutesv', 'cuteSV',
        '--cutesv_params', '"--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL	1000 --diff_ratio_merging_DEL 0.5 --min_support 3 --min_mapq 20 --min_size 30 --max_size -1 --report_readid --genotype "',
        '--output_dir', output_dir,
        '-w', work_dir
    ]
    run_workflow(workflow='long_read_dna_structural_variants.nf',
                 workflow_args=workflow_args)

//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runSqanti3GtfMode } from '../../modules/sqanti3'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.reference_genes_gtf_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf'
params.params_sqanti3_qc = '--report skip'
params.params_sqanti3_filter = 'ml'
params.delete_work_dir = false

if (params.params_sqanti3_qc == true) {
    params_sqanti3_qc = ''
} else {
    params_sqanti3_qc = params.params_sqanti3_qc
}

if (params.params_sqanti3_filter == true) {
    params_sqanti3_filter = ''
} else {
    params_sqanti3_filter = params.params_sqanti3_filter
}

// Step 3. Print inputs and help
log.info """\
         ================================================
         Discover novel isoforms using Sqanti3 (GTF mode)
         ================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run SQANTI3 (GTF mode).

    usage: nexus run --nf-workflow novel_isoform_discovery_sqanti3_gtf_mode.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'gtf_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
        --reference_genes_gtf_file          :   Reference genes GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf).
        --params_sqanti3_qc                 :   sqanti3_qc.py parameters (default: '"--report skip"').
                                                Note that the parameters need to be wrapped in quotes.
        --params_sqanti3_filter             :   sqanti3_filter.py parameters (default: '"ml"').
                                                Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genome_fasta_fai_file     :   ${params.reference_genome_fasta_fai_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_sqanti3_qc                   :   ${params_sqanti3_qc}
        params_sqanti3_filter               :   ${params_sqanti3_filter}
        delete_work_dir                     :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.gtf_file}") }
    .set { input_gtf_files_ch }

// Step 4. Workflow
workflow NOVEL_ISOFORM_DISCOVERY_SQANTI3_GTF_MODE {
    take:
        input_gtf_files_ch
        reference_genome_fasta_file
        reference_genome_fasta_fai_file
        reference_genes_gtf_file
        params_sqanti3_qc
        params_sqanti3_filter
        output_dir
    main:
        runSqanti3GtfMode(
            input_gtf_files_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            reference_genes_gtf_file,
            params_sqanti3_qc,
            params_sqanti3_filter,
            output_dir
        )
    emit:
        runSqanti3GtfMode.out.f
}

workflow {
    NOVEL_ISOFORM_DISCOVERY_SQANTI3_GTF_MODE(
        input_gtf_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_fai_file,
        params.reference_genes_gtf_file,
        params_sqanti3_qc,
        params_sqanti3_filter,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

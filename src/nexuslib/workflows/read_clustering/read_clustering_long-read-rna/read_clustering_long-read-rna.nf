#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { READ_CLUSTERING_GELUSTER }     from '../../../subworkflows/read_clustering/read_clustering_geluster/read_clustering_geluster'
include { READ_CLUSTERING_ISONCLUST3 }   from '../../../subworkflows/read_clustering/read_clustering_isonclust3/read_clustering_isonclust3'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         =========================================================
         Cluster reads in long-read RNA sequencing FASTQ files
         =========================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow read_clustering_long-read-rna.nf -params-file params.yaml [--help]

    All parameters are supplied via a params.yaml file. See params.yaml for
    full documentation and defaults.
    """.stripIndent()
    exit 0
}

// ------------------------------------------------------------
// Step 3. Validate inputs
// ------------------------------------------------------------
def active_methods    = params.methods.toString().tokenize(',').collect { it.trim().toLowerCase() }
def run_all           = active_methods.isEmpty() || active_methods.contains('all')
def run_geluster      = run_all || active_methods.contains('geluster')
def run_isonclust3    = run_all || active_methods.contains('isonclust3')

if (!params.samples_tsv_file) error "ERROR: samples_tsv_file is required."
if (!params.output_dir)       error "ERROR: output_dir is required."

def known_methods = [
    'all',
    'geluster',
    'isonclust3'
]
active_methods.each { m ->
    if (!known_methods.contains(m)) log.warn "WARNING: unknown method '${m}' — will be ignored."
}

log.info """\
    samples_tsv_file             :   ${params.samples_tsv_file}
    output_dir                   :   ${params.output_dir}
    methods                      :   ${params.methods}
    """.stripIndent()

// ------------------------------------------------------------
// Step 4. Set channels
// ------------------------------------------------------------
// FASTQ channel — one gzipped long-read RNA FASTQ per sample
// (used by GeLuster and isONclust3)
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.fastq_file}") }
    .set { input_fastq_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow READ_CLUSTERING_LONGREAD_RNA {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        output_dir
        cfg_geluster
        cfg_isonclust3

    main:
        if (run_geluster) {
            READ_CLUSTERING_GELUSTER(
                input_fastq_files_ch,
                cfg_geluster.extra_args ?: '--rform fq --seqType PacBio',
                output_dir
            )
        }

        if (run_isonclust3) {
            READ_CLUSTERING_ISONCLUST3(
                input_fastq_files_ch,
                cfg_isonclust3.extra_args ?: '--mode pacbio --no-fastq',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    READ_CLUSTERING_LONGREAD_RNA(
        input_fastq_files_ch,
        params.output_dir,
        params.geluster,
        params.isonclust3
    )
}

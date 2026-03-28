#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { decompressFile as decompressGtf }          from '../../../tools/utils'
include { ISOFORM_CHARACTERIZATION_RMATS }           from '../../../subworkflows/isoform_characterization/isoform_characterization_rmats/isoform_characterization_rmats'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ====================================================================
         Characterize isoforms in short-read RNA sequencing BAM files
         ====================================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow isoform_characterization_short-read.nf -params-file params.yaml [--help]

    All parameters are supplied via a params.yaml file. See params.yaml for
    full documentation and defaults.
    """.stripIndent()
    exit 0
}

// ------------------------------------------------------------
// Step 3. Validate inputs
// ------------------------------------------------------------
def active_methods      = params.methods.toString().tokenize(',').collect { it.trim().toLowerCase() }
def run_all             = active_methods.isEmpty() || active_methods.contains('all')
def run_rmats           = run_all || active_methods.contains('rmats')

if (!params.samples_tsv_file)         error "ERROR: samples_tsv_file is required."
if (!params.output_dir)               error "ERROR: output_dir is required."
if (!params.reference_genes_gtf_file) error "ERROR: reference_genes_gtf_file is required."

def known_methods = [
    'all',
    'rmats'
]
active_methods.each { m ->
    if (!known_methods.contains(m)) log.warn "WARNING: unknown method '${m}' — will be ignored."
}

log.info """\
    samples_tsv_file             :   ${params.samples_tsv_file}
    reference_genes_gtf_file     :   ${params.reference_genes_gtf_file}
    output_dir                   :   ${params.output_dir}
    methods                      :   ${params.methods}
    """.stripIndent()

// ------------------------------------------------------------
// Step 4. Set channels
// ------------------------------------------------------------
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow ISOFORM_CHARACTERIZATION_SHORTREAD {
    take:
        input_bam_files_ch          // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genes_gtf_file
        output_dir
        cfg_rmats

    main:
        decompressGtf(reference_genes_gtf_file)
        uncompressed_gtf = decompressGtf.out.f

        if (run_rmats) {
            ISOFORM_CHARACTERIZATION_RMATS(
                input_bam_files_ch,
                uncompressed_gtf,
                cfg_rmats.extra_args ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ISOFORM_CHARACTERIZATION_SHORTREAD(
        input_bam_files_ch,
        params.reference_genes_gtf_file,
        params.output_dir,
        params.rmats
    )
}

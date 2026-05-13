#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { HLA_TYPING_HLAMINER_LR_RNA }   from '../../../subworkflows/hla_typing/hla_typing_hlaminer-lr-rna/hla_typing_hlaminer-lr-rna'
include { HLA_TYPING_SPECIMMUNE }        from '../../../subworkflows/hla_typing/hla_typing_specimmune/hla_typing_specimmune'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ============================================
         HLA typing in long-read sequencing FASTQ
         ============================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow hla_typing_long-read.nf -params-file params.yaml [--help]

    All parameters are supplied via a params.yaml file. See params.yaml for
    full documentation and defaults.
    """.stripIndent()
    exit 0
}

// ------------------------------------------------------------
// Step 3. Validate inputs
// ------------------------------------------------------------
def active_methods       = params.methods.toString().tokenize(',').collect { it.trim().toLowerCase() }
def run_all              = active_methods.isEmpty() || active_methods.contains('all')
def run_hlaminer_lr_rna  = run_all || active_methods.contains('hlaminer-lr-rna')
def run_specimmune       = run_all || active_methods.contains('specimmune')

if (!params.samples_tsv_file) error "ERROR: samples_tsv_file is required."
if (!params.output_dir)       error "ERROR: output_dir is required."

def known_methods = [
    'all',
    'hlaminer-lr-rna',
    'specimmune'
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
// FASTQ channel — single long-read FASTQ per sample (used by both
// HLAminer (lr-rna) and SpecImmune)
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
workflow HLA_TYPING_LONGREAD {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        output_dir
        cfg_hlaminer_lr_rna
        cfg_specimmune

    main:
        if (run_hlaminer_lr_rna) {
            HLA_TYPING_HLAMINER_LR_RNA(
                input_fastq_files_ch,
                cfg_hlaminer_lr_rna.minimap2_extra_args ?: '',
                cfg_hlaminer_lr_rna.hlaminer_extra_args ?: '',
                output_dir
            )
        }

        if (run_specimmune) {
            HLA_TYPING_SPECIMMUNE(
                input_fastq_files_ch,
                cfg_specimmune.extra_args ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    HLA_TYPING_LONGREAD(
        input_fastq_files_ch,
        params.output_dir,
        // YAML key is "hlaminer-lr-rna" (hyphens) — Groovy needs subscript
        // access since the hyphen isn't a valid identifier character.
        params['hlaminer-lr-rna'],
        params.specimmune
    )
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { HLA_TYPING_FUFIHLA }            from '../../../subworkflows/hla_typing/hla_typing_fufihla/hla_typing_fufihla'
include { HLA_TYPING_HLAMINER_LR_DNA }    from '../../../subworkflows/hla_typing/hla_typing_hlaminer-lr-dna/hla_typing_hlaminer-lr-dna'
include { HLA_TYPING_SPECHLA_LONG_READ }  from '../../../subworkflows/hla_typing/hla_typing_spechla-lr-dna/hla_typing_spechla-lr-dna'
include { HLA_TYPING_SPECIMMUNE }         from '../../../subworkflows/hla_typing/hla_typing_specimmune/hla_typing_specimmune'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         =======================================================
         HLA typing in long-read DNA sequencing FASTQ files
         =======================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow hla_typing_long-read-dna.nf -params-file params.yaml [--help]

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
def run_fufihla          = run_all || active_methods.contains('fufihla')
def run_hlaminer_lr_dna  = run_all || active_methods.contains('hlaminer')
def run_spechla_lr_dna   = run_all || active_methods.contains('spechla')
def run_specimmune       = run_all || active_methods.contains('specimmune')

if (!params.samples_tsv_file) error "ERROR: samples_tsv_file is required."
if (!params.output_dir)       error "ERROR: output_dir is required."

def known_methods = [
    'all',
    'fufihla',
    'hlaminer',
    'spechla',
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
// FASTQ channel — single long-read DNA FASTQ per sample
// (used by FuFiHLA, SpecHLA (long-read) and SpecImmune)
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
workflow HLA_TYPING_LONGREAD_DNA {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        output_dir
        cfg_fufihla
        cfg_hlaminer
        cfg_spechla
        cfg_specimmune

    main:
        if (run_fufihla) {
            HLA_TYPING_FUFIHLA(
                input_fastq_files_ch,
                cfg_fufihla.extra_args ?: '--hifi',
                output_dir
            )
        }

        if (run_hlaminer_lr_dna) {
            HLA_TYPING_HLAMINER_LR_DNA(
                input_fastq_files_ch,
                cfg_hlaminer.minimap2_extra_args ?: '-ax map-hifi --secondary=no',
                cfg_hlaminer.hlaminer_extra_args ?: '-s 500 -q 1 -i 1',
                output_dir
            )
        }

        if (run_spechla_lr_dna) {
            HLA_TYPING_SPECHLA_LONG_READ(
                input_fastq_files_ch,
                cfg_spechla.extra_args ?: '',
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
    HLA_TYPING_LONGREAD_DNA(
        input_fastq_files_ch,
        params.output_dir,
        params.fufihla,
        params.hlaminer,
        params.spechla,
        params.specimmune
    )
}

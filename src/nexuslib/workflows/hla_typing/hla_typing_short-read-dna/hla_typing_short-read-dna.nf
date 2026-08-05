#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { HLA_TYPING_HLAMINER_SR_DNA }     from '../../../subworkflows/hla_typing/hla_typing_hlaminer-sr-dna/hla_typing_hlaminer-sr-dna'
include { HLA_TYPING_OPTITYPE }            from '../../../subworkflows/hla_typing/hla_typing_optitype/hla_typing_optitype'
include { HLA_TYPING_SPECHLA_SHORT_READ }  from '../../../subworkflows/hla_typing/hla_typing_spechla/hla_typing_spechla'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         =======================================================
         HLA typing in short-read DNA sequencing FASTQ files
         =======================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow hla_typing_short-read-dna.nf -params-file params.yaml [--help]

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
def run_hlaminer_sr_dna = run_all || active_methods.contains('hlaminer')
def run_optitype        = run_all || active_methods.contains('optitype')
def run_spechla         = run_all || active_methods.contains('spechla')

if (!params.samples_tsv_file) error "ERROR: samples_tsv_file is required."
if (!params.output_dir)       error "ERROR: output_dir is required."

def known_methods = [
    'all',
    'hlaminer',
    'optitype',
    'spechla'
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
// FASTQ channel — paired-end short-read DNA (used by OptiType and SpecHLA)
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.fastq_file_1}",
        "${row.fastq_file_2}") }
    .set { input_fastq_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow HLA_TYPING_SHORTREAD_DNA {
    take:
        input_fastq_files_ch        // channel: [val(sample_id), path(fastq_file_1), path(fastq_file_2)]
        output_dir
        cfg_hlaminer
        cfg_optitype
        cfg_spechla

    main:
        if (run_hlaminer_sr_dna) {
            HLA_TYPING_HLAMINER_SR_DNA(
                input_fastq_files_ch,
                cfg_hlaminer.bwa_aln_extra_args ?: '-e 0 -o 0',
                cfg_hlaminer.bwa_sampe_extra_args ?: '-o 1000',
                cfg_hlaminer.hlaminer_extra_args ?: '-s 500',
                output_dir
            )
        }

        if (run_optitype) {
            HLA_TYPING_OPTITYPE(
                input_fastq_files_ch,
                cfg_optitype.extra_args ?: '--dna',
                output_dir
            )
        }

        if (run_spechla) {
            HLA_TYPING_SPECHLA_SHORT_READ(
                input_fastq_files_ch,
                cfg_spechla.extra_args ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    HLA_TYPING_SHORTREAD_DNA(
        input_fastq_files_ch,
        params.output_dir,
        params.hlaminer,
        params.optitype,
        params.spechla
    )
}

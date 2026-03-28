#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { HLA_TYPING_ARCASHLA }     from '../../../subworkflows/hla_typing/hla_typing_arcashla/hla_typing_arcashla'
include { HLA_TYPING_HLAPROFILER }  from '../../../subworkflows/hla_typing/hla_typing_hlaprofiler/hla_typing_hlaprofiler'
include { HLA_TYPING_SEQ2HLA }      from '../../../subworkflows/hla_typing/hla_typing_seq2hla/hla_typing_seq2hla'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         =======================================================
         HLA typing in short-read RNA sequencing BAM/FASTQ files
         =======================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow hla_typing_short-read.nf -params-file params.yaml [--help]

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
def run_arcashla        = run_all || active_methods.contains('arcashla')
def run_hlaprofiler     = run_all || active_methods.contains('hlaprofiler')
def run_seq2hla         = run_all || active_methods.contains('seq2hla')

if (!params.samples_tsv_file) error "ERROR: samples_tsv_file is required."
if (!params.output_dir)       error "ERROR: output_dir is required."

def known_methods = [
    'all',
    'arcashla',
    'hlaprofiler',
    'seq2hla'
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
// BAM channel for arcasHLA
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// FASTQ channel for HLAProfiler and seq2HLA
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
workflow HLA_TYPING_SHORTREAD {
    take:
        input_bam_files_ch          // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        input_fastq_files_ch        // channel: [val(sample_id), path(fastq_file_1), path(fastq_file_2)]
        output_dir
        cfg_hlaprofiler
        cfg_seq2hla

    main:
        if (run_arcashla) {
            HLA_TYPING_ARCASHLA(
                input_bam_files_ch,
                output_dir
            )
        }

        if (run_hlaprofiler) {
            HLA_TYPING_HLAPROFILER(
                input_fastq_files_ch,
                cfg_hlaprofiler.extra_args ?: '',
                output_dir
            )
        }

        if (run_seq2hla) {
            HLA_TYPING_SEQ2HLA(
                input_fastq_files_ch,
                cfg_seq2hla.extra_args ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    HLA_TYPING_SHORTREAD(
        input_bam_files_ch,
        input_fastq_files_ch,
        params.output_dir,
        params.hlaprofiler,
        params.seq2hla
    )
}

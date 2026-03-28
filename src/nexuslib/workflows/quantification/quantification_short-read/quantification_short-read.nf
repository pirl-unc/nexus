#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { QUANTIFICATION_KALLISTO_PE }              from '../../../subworkflows/quantification/quantification_kallisto-pe/quantification_kallisto-pe'
include { QUANTIFICATION_SALMON_FASTQ }             from '../../../subworkflows/quantification/quantification_salmon-fastq/quantification_salmon-fastq'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ================================================================
         Quantify RNA in short-read RNA sequencing FASTQ files
         ================================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow quantification_short-read.nf -params-file params.yaml [--help]

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
def run_kallisto        = run_all || active_methods.contains('kallisto')
def run_salmon          = run_all || active_methods.contains('salmon')

if (!params.samples_tsv_file)                error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                      error "ERROR: output_dir is required."
if (!params.reference_transcripts_fasta_file) error "ERROR: reference_transcripts_fasta_file is required."

def known_methods = [
    'all',
    'kallisto',
    'salmon'
]
active_methods.each { m ->
    if (!known_methods.contains(m)) log.warn "WARNING: unknown method '${m}' — will be ignored."
}

log.info """\
    samples_tsv_file                 :   ${params.samples_tsv_file}
    reference_transcripts_fasta_file :   ${params.reference_transcripts_fasta_file}
    output_dir                       :   ${params.output_dir}
    methods                          :   ${params.methods}
    """.stripIndent()

// ------------------------------------------------------------
// Step 4. Set channels
// ------------------------------------------------------------
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
workflow QUANTIFICATION_SHORTREAD {
    take:
        input_fastq_files_ch        // channel: [val(sample_id), path(fastq_file_1), path(fastq_file_2)]
        reference_transcripts_fasta_file
        output_dir
        cfg_kallisto
        cfg_salmon

    main:
        if (run_kallisto) {
            QUANTIFICATION_KALLISTO_PE(
                input_fastq_files_ch,
                reference_transcripts_fasta_file,
                cfg_kallisto.quant_extra_args ?: '',
                output_dir
            )
        }

        if (run_salmon) {
            QUANTIFICATION_SALMON_FASTQ(
                input_fastq_files_ch,
                reference_transcripts_fasta_file,
                cfg_salmon.index_extra_args ?: '',
                cfg_salmon.quant_extra_args ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    QUANTIFICATION_SHORTREAD(
        input_fastq_files_ch,
        params.reference_transcripts_fasta_file,
        params.output_dir,
        params.kallisto,
        params.salmon
    )
}

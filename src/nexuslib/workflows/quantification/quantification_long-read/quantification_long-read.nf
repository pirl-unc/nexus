#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { decompressFile as decompressFasta }       from '../../../tools/utils'
include { decompressFile as decompressGtf }         from '../../../tools/utils'
include { decompressFile as decompressTranscripts } from '../../../tools/utils'
include { QUANTIFICATION_BAMBU }                    from '../../../subworkflows/quantification/quantification_bambu/quantification_bambu'
include { QUANTIFICATION_KALLISTO_LR }              from '../../../subworkflows/quantification/quantification_kallisto-lr/quantification_kallisto-lr'
include { QUANTIFICATION_LIQA }                     from '../../../subworkflows/quantification/quantification_liqa/quantification_liqa'
include { QUANTIFICATION_OARFISH }                  from '../../../subworkflows/quantification/quantification_oarfish/quantification_oarfish'
include { QUANTIFICATION_TRANSIGNER }               from '../../../subworkflows/quantification/quantification_transigner/quantification_transigner'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ================================================================
         Quantify RNA in long-read RNA sequencing FASTQ/BAM files
         ================================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow quantification_long-read.nf -params-file params.yaml [--help]

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
def run_bambu           = run_all || active_methods.contains('bambu')
def run_kallisto        = run_all || active_methods.contains('kallisto')
def run_liqa            = run_all || active_methods.contains('liqa')
def run_oarfish         = run_all || active_methods.contains('oarfish')
def run_transigner      = run_all || active_methods.contains('transigner')

if (!params.samples_tsv_file)            error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                  error "ERROR: output_dir is required."

if (run_bambu) {
    if (!params.reference_genome_fasta_file) error "ERROR: reference_genome_fasta_file is required when running bambu."
    if (!params.reference_genes_gtf_file)    error "ERROR: reference_genes_gtf_file is required when running bambu."
}

if (run_kallisto) {
    if (!params.reference_transcripts_fasta_file) error "ERROR: reference_transcripts_fasta_file is required when running kallisto."
    if (!params.reference_genes_gtf_file)         error "ERROR: reference_genes_gtf_file is required when running kallisto."
}

if (run_liqa) {
    if (!params.reference_genes_gtf_file) error "ERROR: reference_genes_gtf_file is required when running liqa."
}

if (run_oarfish) {
    if (!params.reference_transcripts_fasta_file) error "ERROR: reference_transcripts_fasta_file is required when running oarfish."
}

if (run_transigner) {
    if (!params.reference_transcripts_fasta_file) error "ERROR: reference_transcripts_fasta_file is required when running transigner."
}

def known_methods = [
    'all',
    'bambu',
    'kallisto',
    'liqa',
    'oarfish',
    'transigner'
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
// BAM channel for bambu and liqa
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// FASTQ channel for kallisto, oarfish, transigner
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
workflow QUANTIFICATION_LONGREAD {
    take:
        input_bam_files_ch          // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        input_fastq_files_ch        // channel: [val(sample_id), path(fastq_file)]
        reference_genome_fasta_file
        reference_genes_gtf_file
        reference_transcripts_fasta_file
        output_dir
        cfg_bambu
        cfg_kallisto
        cfg_liqa
        cfg_oarfish
        cfg_transigner

    main:
        if (run_bambu) {
            QUANTIFICATION_BAMBU(
                input_bam_files_ch,
                reference_genome_fasta_file,
                reference_genes_gtf_file,
                output_dir
            )
        }

        if (run_kallisto) {
            QUANTIFICATION_KALLISTO_LR(
                input_fastq_files_ch,
                reference_transcripts_fasta_file,
                reference_genes_gtf_file,
                cfg_kallisto.index_extra_args ?: '',
                cfg_kallisto.bus_extra_args ?: '',
                cfg_kallisto.bustools_sort_extra_args ?: '',
                cfg_kallisto.bustools_count_extra_args ?: '',
                cfg_kallisto.quanttcc_extra_args ?: '',
                output_dir
            )
        }

        if (run_liqa) {
            QUANTIFICATION_LIQA(
                input_bam_files_ch,
                reference_genes_gtf_file,
                cfg_liqa.extra_args ?: '',
                output_dir
            )
        }

        if (run_oarfish) {
            QUANTIFICATION_OARFISH(
                input_fastq_files_ch,
                reference_transcripts_fasta_file,
                cfg_oarfish.extra_args ?: '',
                output_dir
            )
        }

        if (run_transigner) {
            QUANTIFICATION_TRANSIGNER(
                input_fastq_files_ch,
                reference_transcripts_fasta_file,
                cfg_transigner.align_extra_args ?: '',
                cfg_transigner.prefilter_extra_args ?: '',
                cfg_transigner.em_extra_args ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    QUANTIFICATION_LONGREAD(
        input_bam_files_ch,
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params.reference_transcripts_fasta_file,
        params.output_dir,
        params.bambu,
        params.kallisto,
        params.liqa,
        params.oarfish,
        params.transigner
    )
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { VARIANT_CALLING_CLAIR3RNA }    from '../../../subworkflows/variant_calling/variant_calling_clair3rna/variant_calling_clair3rna'
include { VARIANT_CALLING_DE_SOUZA }     from '../../../subworkflows/variant_calling/variant_calling_de-souza/variant_calling_de-souza'
include { VARIANT_CALLING_LONGGF }       from '../../../subworkflows/variant_calling/variant_calling_longgf/variant_calling_longgf'
include { VARIANT_CALLING_PBFUSION }     from '../../../subworkflows/variant_calling/variant_calling_pbfusion/variant_calling_pbfusion'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ===============================================================
         Identify RNA variants and gene fusions in long-read RNA samples
         ===============================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow variant_calling_long-read-rna.nf -params-file params.yaml [--help]

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
def run_clair3rna       = run_all || active_methods.contains('clair3rna')
def run_de_souza        = run_all || active_methods.contains('de-souza') || active_methods.contains('de_souza')
def run_longgf          = run_all || active_methods.contains('longgf')
def run_pbfusion        = run_all || active_methods.contains('pbfusion')

if (!params.samples_tsv_file)            error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                  error "ERROR: output_dir is required."
if (!params.reference_genome_fasta_file) error "ERROR: reference_genome_fasta_file is required."

if (run_longgf || run_pbfusion) {
    if (!params.reference_genes_gtf_file) error "ERROR: reference_genes_gtf_file is required when running longgf or pbfusion."
}

if (run_de_souza) {
    if (!params.de_souza.deepvariant_input_path)  error "ERROR: de_souza.deepvariant_input_path is required when running de-souza."
    if (!params.de_souza.deepvariant_output_path) error "ERROR: de_souza.deepvariant_output_path is required when running de-souza."
}

def known_methods = [
    'all',
    'clair3rna',
    'de-souza',
    'longgf',
    'pbfusion'
]
active_methods.each { m ->
    if (!known_methods.contains(m)) log.warn "WARNING: unknown method '${m}' — will be ignored."
}

log.info """\
    samples_tsv_file             :   ${params.samples_tsv_file}
    reference_genome_fasta_file  :   ${params.reference_genome_fasta_file}
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
        "${row.bam_bai_file}",
        "${row.fastq_file}") }
    .set { input_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_LONGREAD_RNA {
    take:
        input_files_ch              // channel: [val(sample_id), path(bam_file), path(bam_bai_file), path(fastq_file)]
        reference_genome_fasta_file
        reference_genes_gtf_file
        output_dir
        cfg_clair3rna
        cfg_de_souza
        cfg_longgf
        cfg_pbfusion

    main:
        // clair3rna, longgf, and pbfusion take an aligned BAM;
        // de-souza takes a raw FASTQ (it aligns internally with minimap2).
        bam_files_ch = input_files_ch.map { it ->
            tuple(it[0], it[1], it[2])
        }
        fastq_files_ch = input_files_ch.map { it ->
            tuple(it[0], it[3])
        }

        if (run_clair3rna) {
            VARIANT_CALLING_CLAIR3RNA(
                bam_files_ch,
                reference_genome_fasta_file,
                cfg_clair3rna.extra_args ?: '',
                output_dir
            )
        }

        if (run_de_souza) {
            VARIANT_CALLING_DE_SOUZA(
                fastq_files_ch,
                reference_genome_fasta_file,
                cfg_de_souza.params_minimap2 ?: '',
                cfg_de_souza.platform_tag ?: 'unknown',
                cfg_de_souza.platform_unit_tag ?: 'unknown',
                cfg_de_souza.library_tag ?: 'unknown',
                cfg_de_souza.params_samtools_view ?: '',
                cfg_de_souza.deepvariant_containerization,
                cfg_de_souza.deepvariant_bin_version,
                cfg_de_souza.deepvariant_bin_path,
                cfg_de_souza.deepvariant_input_path,
                cfg_de_souza.deepvariant_output_path,
                cfg_de_souza.deepvariant_model_type,
                output_dir
            )
        }

        if (run_longgf) {
            VARIANT_CALLING_LONGGF(
                bam_files_ch,
                reference_genes_gtf_file,
                cfg_longgf.extra_args ?: '',
                output_dir
            )
        }

        if (run_pbfusion) {
            VARIANT_CALLING_PBFUSION(
                bam_files_ch,
                reference_genome_fasta_file,
                reference_genes_gtf_file,
                cfg_pbfusion.extra_args ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_LONGREAD_RNA(
        input_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params.output_dir,
        params.clair3rna,
        params.de_souza,
        params.longgf,
        params.pbfusion
    )
}

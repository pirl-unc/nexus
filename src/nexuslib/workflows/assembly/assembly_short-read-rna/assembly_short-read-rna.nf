#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { ASSEMBLY_SPADES }     from '../../../subworkflows/assembly/assembly_spades/assembly_spades'
include { ASSEMBLY_TRINITY }    from '../../../subworkflows/assembly/assembly_trinity/assembly_trinity'
include { ASSEMBLY_BOOKEND }    from '../../../subworkflows/assembly/assembly_bookend/assembly_bookend'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ==================================================
         Assemble transcripts from short-read RNA samples
         ==================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow assembly_short-read-rna.nf -params-file params.yaml [--help]

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
def run_spades          = run_all || active_methods.contains('spades')
def run_trinity         = run_all || active_methods.contains('trinity')
def run_bookend         = run_all || active_methods.contains('bookend')

if (!params.samples_tsv_file)            error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                  error "ERROR: output_dir is required."

if (run_bookend) {
    if (!params.reference_genome_fasta_file) error "ERROR: reference_genome_fasta_file is required when running bookend."
}

def known_methods = [
    'all',
    'spades',
    'trinity',
    'bookend'
]
active_methods.each { m ->
    if (!known_methods.contains(m)) log.warn "WARNING: unknown method '${m}' — will be ignored."
}

log.info """\
    samples_tsv_file             :   ${params.samples_tsv_file}
    reference_genome_fasta_file  :   ${params.reference_genome_fasta_file}
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
        "${row.fastq_file_1}",
        "${row.fastq_file_2}",
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow ASSEMBLY_SHORTREAD_RNA {
    take:
        input_files_ch              // channel: [val(sample_id), path(fastq_file_1), path(fastq_file_2), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        output_dir
        cfg_spades
        cfg_trinity
        cfg_bookend

    main:
        // spades and trinity assemble de novo from paired-end FASTQ;
        // bookend assembles from an aligned BAM (plus the reference FASTA).
        fastq_files_ch = input_files_ch.map { it ->
            tuple(it[0], it[1], it[2])
        }
        bam_files_ch = input_files_ch.map { it ->
            tuple(it[0], it[3], it[4])
        }

        if (run_spades) {
            ASSEMBLY_SPADES(
                fastq_files_ch,
                cfg_spades.extra_args ?: '',
                output_dir
            )
        }

        if (run_trinity) {
            ASSEMBLY_TRINITY(
                fastq_files_ch,
                cfg_trinity.extra_args ?: '',
                output_dir
            )
        }

        if (run_bookend) {
            ASSEMBLY_BOOKEND(
                bam_files_ch,
                reference_genome_fasta_file,
                cfg_bookend.params_assemble ?: '',
                cfg_bookend.params_fasta ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ASSEMBLY_SHORTREAD_RNA(
        input_files_ch,
        params.reference_genome_fasta_file,
        params.output_dir,
        params.spades,
        params.trinity,
        params.bookend
    )
}

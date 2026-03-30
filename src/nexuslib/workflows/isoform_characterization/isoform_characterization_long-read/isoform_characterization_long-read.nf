#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { decompressFile as decompressFasta }        from '../../../tools/utils'
include { decompressFile as decompressGtf }          from '../../../tools/utils'
include { extractGtfFromDir as extractEspressoGtf }    from '../../../tools/utils'
include { extractGtfFromDir as extractIsoquantGtf }    from '../../../tools/utils'
include { extractGtfFromDir as extractIsotoolsGtf }    from '../../../tools/utils'
include { extractGtfFromDir as extractMandalorionGtf } from '../../../tools/utils'
include { ISOFORM_CHARACTERIZATION_ESPRESSO }        from '../../../subworkflows/isoform_characterization/isoform_characterization_espresso/isoform_characterization_espresso'
include { ISOFORM_CHARACTERIZATION_FLAIR }           from '../../../subworkflows/isoform_characterization/isoform_characterization_flair/isoform_characterization_flair'
include { ISOFORM_CHARACTERIZATION_ISOQUANT }        from '../../../subworkflows/isoform_characterization/isoform_characterization_isoquant/isoform_characterization_isoquant'
include { ISOFORM_CHARACTERIZATION_ISOTOOLS }        from '../../../subworkflows/isoform_characterization/isoform_characterization_isotools/isoform_characterization_isotools'
include { ISOFORM_CHARACTERIZATION_MANDALORION }     from '../../../subworkflows/isoform_characterization/isoform_characterization_mandalorion/isoform_characterization_mandalorion'
include { ISOFORM_CHARACTERIZATION_SQANTI3_GTF as SQANTI3_GTF_ESPRESSO }     from '../../../subworkflows/isoform_characterization/isoform_characterization_sqanti3-gtf/isoform_characterization_sqanti3-gtf'
include { ISOFORM_CHARACTERIZATION_SQANTI3_GTF as SQANTI3_GTF_FLAIR }        from '../../../subworkflows/isoform_characterization/isoform_characterization_sqanti3-gtf/isoform_characterization_sqanti3-gtf'
include { ISOFORM_CHARACTERIZATION_SQANTI3_GTF as SQANTI3_GTF_ISOQUANT }     from '../../../subworkflows/isoform_characterization/isoform_characterization_sqanti3-gtf/isoform_characterization_sqanti3-gtf'
include { ISOFORM_CHARACTERIZATION_SQANTI3_GTF as SQANTI3_GTF_ISOTOOLS }     from '../../../subworkflows/isoform_characterization/isoform_characterization_sqanti3-gtf/isoform_characterization_sqanti3-gtf'
include { ISOFORM_CHARACTERIZATION_SQANTI3_GTF as SQANTI3_GTF_MANDALORION }  from '../../../subworkflows/isoform_characterization/isoform_characterization_sqanti3-gtf/isoform_characterization_sqanti3-gtf'
include { ISOFORM_CHARACTERIZATION_TALON }           from '../../../subworkflows/isoform_characterization/isoform_characterization_talon/isoform_characterization_talon'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         =================================================================
         Characterize isoforms in long-read RNA sequencing BAM/FASTQ files
         =================================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow isoform_characterization_long-read.nf -params-file params.yaml [--help]

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
def run_espresso        = run_all || active_methods.contains('espresso')
def run_flair           = run_all || active_methods.contains('flair')
def run_isoquant        = run_all || active_methods.contains('isoquant')
def run_isotools        = run_all || active_methods.contains('isotools')
def run_mandalorion     = run_all || active_methods.contains('mandalorion')
def run_sqanti3_gtf     = run_all || active_methods.contains('sqanti3-gtf')
def run_talon           = run_all || active_methods.contains('talon')

if (!params.samples_tsv_file)            error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                  error "ERROR: output_dir is required."
if (!params.reference_genome_fasta_file) error "ERROR: reference_genome_fasta_file is required."
if (!params.reference_genes_gtf_file)    error "ERROR: reference_genes_gtf_file is required."

def known_methods = [
    'all',
    'espresso',
    'flair',
    'isoquant',
    'isotools',
    'mandalorion',
    'sqanti3-gtf',
    'talon'
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
// BAM channel for ESPRESSO, IsoTools, TALON
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// FASTQ channel for FLAIR, IsoQuant, Mandalorion
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
workflow ISOFORM_CHARACTERIZATION_LONGREAD {
    take:
        input_bam_files_ch          // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        input_fastq_files_ch        // channel: [val(sample_id), path(fastq_file)]
        reference_genome_fasta_file
        reference_genes_gtf_file
        output_dir
        cfg_espresso
        cfg_flair
        cfg_isoquant
        cfg_isotools
        cfg_mandalorion
        cfg_sqanti3_gtf
        cfg_talon

    main:
        decompressFasta(reference_genome_fasta_file)
        decompressGtf(reference_genes_gtf_file)
        uncompressed_fasta = decompressFasta.out.f
        uncompressed_gtf = decompressGtf.out.f

        if (run_espresso) {
            ISOFORM_CHARACTERIZATION_ESPRESSO(
                input_bam_files_ch,
                uncompressed_fasta,
                uncompressed_gtf,
                cfg_espresso.s_extra_args ?: '',
                cfg_espresso.c_extra_args ?: '',
                cfg_espresso.q_extra_args ?: '',
                output_dir
            )
            if (run_sqanti3_gtf) {
                extractEspressoGtf(ISOFORM_CHARACTERIZATION_ESPRESSO.out, "*_updated.gtf")
                SQANTI3_GTF_ESPRESSO(
                    extractEspressoGtf.out.f,
                    uncompressed_fasta,
                    uncompressed_gtf,
                    cfg_sqanti3_gtf.qc_extra_args ?: '',
                    cfg_sqanti3_gtf.filter_extra_args ?: '',
                    "espresso",
                    output_dir
                )
            }
        }

        if (run_flair) {
            ISOFORM_CHARACTERIZATION_FLAIR(
                input_fastq_files_ch,
                uncompressed_fasta,
                uncompressed_gtf,
                cfg_flair.align_extra_args ?: '',
                cfg_flair.correct_extra_args ?: '',
                cfg_flair.collapse_extra_args ?: '',
                output_dir
            )
            if (run_sqanti3_gtf) {
                flair_gtf_ch = ISOFORM_CHARACTERIZATION_FLAIR.out
                    .map { sid, bed, fasta, gtf -> tuple(sid, gtf) }
                SQANTI3_GTF_FLAIR(
                    flair_gtf_ch,
                    uncompressed_fasta,
                    uncompressed_gtf,
                    cfg_sqanti3_gtf.qc_extra_args ?: '',
                    cfg_sqanti3_gtf.filter_extra_args ?: '',
                    "flair",
                    output_dir
                )
            }
        }

        if (run_isoquant) {
            ISOFORM_CHARACTERIZATION_ISOQUANT(
                input_fastq_files_ch,
                uncompressed_fasta,
                uncompressed_gtf,
                cfg_isoquant.extra_args ?: '',
                output_dir
            )
            if (run_sqanti3_gtf) {
                extractIsoquantGtf(ISOFORM_CHARACTERIZATION_ISOQUANT.out, "*.transcript_models.gtf")
                SQANTI3_GTF_ISOQUANT(
                    extractIsoquantGtf.out.f,
                    uncompressed_fasta,
                    uncompressed_gtf,
                    cfg_sqanti3_gtf.qc_extra_args ?: '',
                    cfg_sqanti3_gtf.filter_extra_args ?: '',
                    "isoquant",
                    output_dir
                )
            }
        }

        if (run_isotools) {
            ISOFORM_CHARACTERIZATION_ISOTOOLS(
                input_bam_files_ch,
                uncompressed_fasta,
                uncompressed_gtf,
                cfg_isotools.extra_args ?: '',
                output_dir
            )
            if (run_sqanti3_gtf) {
                extractIsotoolsGtf(ISOFORM_CHARACTERIZATION_ISOTOOLS.out, "*.gtf")
                SQANTI3_GTF_ISOTOOLS(
                    extractIsotoolsGtf.out.f,
                    uncompressed_fasta,
                    uncompressed_gtf,
                    cfg_sqanti3_gtf.qc_extra_args ?: '',
                    cfg_sqanti3_gtf.filter_extra_args ?: '',
                    "isotools",
                    output_dir
                )
            }
        }

        if (run_mandalorion) {
            ISOFORM_CHARACTERIZATION_MANDALORION(
                input_fastq_files_ch,
                uncompressed_fasta,
                uncompressed_gtf,
                cfg_mandalorion.extra_args ?: '',
                output_dir
            )
            if (run_sqanti3_gtf) {
                extractMandalorionGtf(ISOFORM_CHARACTERIZATION_MANDALORION.out, "*.gtf")
                SQANTI3_GTF_MANDALORION(
                    extractMandalorionGtf.out.f,
                    uncompressed_fasta,
                    uncompressed_gtf,
                    cfg_sqanti3_gtf.qc_extra_args ?: '',
                    cfg_sqanti3_gtf.filter_extra_args ?: '',
                    "mandalorion",
                    output_dir
                )
            }
        }

        if (run_talon) {
            ISOFORM_CHARACTERIZATION_TALON(
                input_bam_files_ch,
                uncompressed_gtf,
                cfg_talon.initdb_extra_args ?: '',
                cfg_talon.extra_args ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ISOFORM_CHARACTERIZATION_LONGREAD(
        input_bam_files_ch,
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params.output_dir,
        params.espresso,
        params.flair,
        params.isoquant,
        params.isotools,
        params.mandalorion,
        params.sqanti3_gtf,
        params.talon
    )
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { HAPLOTYPE_PHASING_FLAIR_LONGSHOT }    from '../../../subworkflows/haplotype_phasing/haplotype_phasing_flair-longshot/haplotype_phasing_flair-longshot'
include { HAPLOTYPE_PHASING_LONGCALLR }         from '../../../subworkflows/haplotype_phasing/haplotype_phasing_longcallr/haplotype_phasing_longcallr'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ===============================================================================
         Phase variants in long-read RNA sequencing data
         (FLAIR + Longshot, LongcallR)
         ===============================================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow haplotype_phasing_long-read-rna.nf -params-file params.yaml [--help]

    All parameters are supplied via a params.yaml file. See params.yaml for
    full documentation and defaults.
    """.stripIndent()
    exit 0
}

// ------------------------------------------------------------
// Step 3. Validate inputs
// ------------------------------------------------------------
def active_methods   = params.methods.toString().tokenize(',').collect { it.trim().toLowerCase() }
def run_all          = active_methods.isEmpty() || active_methods.contains('all')
def run_flair_ls     = run_all || active_methods.contains('flair-longshot')
def run_longcallr    = run_all || active_methods.contains('longcallr')

if (!params.samples_tsv_file)            error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                  error "ERROR: output_dir is required."
if (!params.reference_genome_fasta_file) error "ERROR: reference_genome_fasta_file is required."

if (run_longcallr) {
    if (!params.longcallr.preset)
        error "ERROR: longcallr.preset is required when running longcallr (choices: hifi-isoseq, hifi-masseq, ont-cdna, ont-drna)."
    if (!params.longcallr.reference_genes_gtf_file)
        error "ERROR: longcallr.reference_genes_gtf_file is required when running longcallr."

    def known_presets = ['hifi-isoseq', 'hifi-masseq', 'ont-cdna', 'ont-drna']
    if (!known_presets.contains(params.longcallr.preset.toString().toLowerCase()))
        error "ERROR: longcallr.preset must be one of ${known_presets} (got: '${params.longcallr.preset}')."
}

def known_methods = ['all', 'flair-longshot', 'longcallr']
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
// FASTQ channel — required when running flair-longshot.
//   samples.tsv columns: sample_id, fastq_file
if (run_flair_ls) {
    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_fastq_files_ch }
}

// BAM channel — required when running longcallr.
//   samples.tsv columns: sample_id, bam_file, bam_bai_file
if (run_longcallr) {
    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }
}

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow HAPLOTYPE_PHASING_LONGREAD_RNA {
    take:
        input_fastq_files_ch        // channel: [val(sample_id), path(fastq_file)] — used by flair-longshot
        input_bam_files_ch          // channel: [val(sample_id), path(bam_file), path(bam_bai_file)] — used by longcallr
        reference_genome_fasta_file
        output_dir
        cfg_flair_longshot
        cfg_longcallr

    main:
        if (run_flair_ls) {
            HAPLOTYPE_PHASING_FLAIR_LONGSHOT(
                input_fastq_files_ch,
                reference_genome_fasta_file,
                cfg_flair_longshot.flair_align_extra_args ?: '',
                cfg_flair_longshot.longshot_extra_args    ?: '',
                output_dir
            )
        }

        if (run_longcallr) {
            HAPLOTYPE_PHASING_LONGCALLR(
                input_bam_files_ch,
                cfg_longcallr.preset,
                reference_genome_fasta_file,
                cfg_longcallr.reference_genes_gtf_file,
                cfg_longcallr.extra_args ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    HAPLOTYPE_PHASING_LONGREAD_RNA(
        run_flair_ls   ? input_fastq_files_ch : Channel.empty(),
        run_longcallr  ? input_bam_files_ch   : Channel.empty(),
        params.reference_genome_fasta_file,
        params.output_dir,
        params.flair_longshot,
        params.longcallr
    )
}

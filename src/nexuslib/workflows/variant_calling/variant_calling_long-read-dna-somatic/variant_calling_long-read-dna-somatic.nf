#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { decompressFile as decompressFasta }           from '../../../tools/utils'
include { VARIANT_CALLING_CLAIRS }                      from '../../../subworkflows/variant_calling/variant_calling_clairs/variant_calling_clairs'
include { VARIANT_CALLING_DEEPSOMATIC }                 from '../../../subworkflows/variant_calling/variant_calling_deepsomatic/variant_calling_deepsomatic'
include { VARIANT_CALLING_DELLY2_LONGREAD_SOMATIC }     from '../../../subworkflows/variant_calling/variant_calling_delly2-lr-somatic/variant_calling_delly2-lr-somatic'
include { VARIANT_CALLING_DYSGU_SOMATIC }               from '../../../subworkflows/variant_calling/variant_calling_dysgu-somatic/variant_calling_dysgu-somatic'
include { VARIANT_CALLING_NANOMONSV }                   from '../../../subworkflows/variant_calling/variant_calling_nanomonsv/variant_calling_nanomonsv'
include { VARIANT_CALLING_SAVANA }                      from '../../../subworkflows/variant_calling/variant_calling_savana/variant_calling_savana'
include { VARIANT_CALLING_SEVERUS }                     from '../../../subworkflows/variant_calling/variant_calling_severus/variant_calling_severus'
include { VARIANT_PHASING_WHATSHAP }                    from '../../../subworkflows/variant_phasing/variant_phasing_whatshap/variant_phasing_whatshap'
include { VARIANT_CALLING_SVISIONPRO }                  from '../../../subworkflows/variant_calling/variant_calling_svisionpro/variant_calling_svisionpro'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ===================================================================
         Identify somatic DNA variants in long-read DNA sequencing BAM files
         ===================================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow variant_calling_long-read-dna-somatic.nf -params-file params.yaml [--help]

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
def run_clairs          = run_all || active_methods.contains('clairs')
def run_deepsomatic     = run_all || active_methods.contains('deepsomatic')
def run_delly2          = run_all || active_methods.contains('delly2')
def run_dysgu           = run_all || active_methods.contains('dysgu')
def run_nanomonsv       = run_all || active_methods.contains('nanomonsv')
def run_savana          = run_all || active_methods.contains('savana')
def run_severus         = run_all || active_methods.contains('severus')
def run_svisionpro      = run_all || active_methods.contains('svisionpro')

if (!params.samples_tsv_file)            error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                  error "ERROR: output_dir is required."
if (!params.reference_genome_fasta_file) error "ERROR: reference_genome_fasta_file is required."

if (run_deepsomatic) {
    if (!params.deepsomatic.input_path)  error "ERROR: deepsomatic.input_path is required when running deepsomatic."
    if (!params.deepsomatic.output_path) error "ERROR: deepsomatic.output_path is required when running deepsomatic."
}

if (run_delly2) {
    if (!params.delly2.exclude_tsv_file) error "ERROR: delly2.exclude_tsv_file is required when running delly2."
}

if (run_savana) {
    if (!params.savana.contigs_txt_file)   error "ERROR: savana.contigs_txt_file is required when running savana."
    if (!params.savana.custom_params_file) error "ERROR: savana.custom_params_file is required when running savana."
}

if (run_severus) {
    if (!params.severus.vntr_bed_file) error "ERROR: severus.vntr_bed_file is required when running severus."
}

if (run_svisionpro) {
    if (!params.svisionpro.model_file) error "ERROR: svisionpro.model_file is required when running svisionpro."
}

def known_methods = [
    'all',
    'clairs',
    'deepsomatic',
    'delly2',
    'dysgu',
    'nanomonsv',
    'savana',
    'severus',
    'svisionpro'
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
// Standard 5-element channel for most callers
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.tumor_bam_file}",
        "${row.tumor_bam_bai_file}",
        "${row.normal_bam_file}",
        "${row.normal_bam_bai_file}") }
    .set { input_bam_files_ch }

// Channel for WhatsHap phasing of normal small variant VCF (for Severus)
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.normal_bam_file}",
        "${row.normal_bam_bai_file}",
        "${row.normal_small_variants_vcf_file}") }
    .set { whatshap_input_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_LONGREAD_SOMATIC {
    take:
        input_bam_files_ch          // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)]
        whatshap_input_ch           // channel: [val(sample_id), path(normal_bam_file), path(normal_bam_bai_file), path(normal_small_variants_vcf_file)] — for WhatsHap → Severus
        reference_genome_fasta_file
        output_dir
        cfg_clairs
        cfg_deepsomatic
        cfg_delly2
        cfg_dysgu
        cfg_nanomonsv
        cfg_savana
        cfg_severus
        cfg_svisionpro

    main:
        decompressFasta(reference_genome_fasta_file)
        uncompressed_fasta = decompressFasta.out.f

        if (run_clairs) {
            VARIANT_CALLING_CLAIRS(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_clairs.extra_args ?: '',
                output_dir
            )
        }

        if (run_deepsomatic) {
            VARIANT_CALLING_DEEPSOMATIC(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_deepsomatic.containerization,
                cfg_deepsomatic.bin_version,
                cfg_deepsomatic.bin_path,
                cfg_deepsomatic.input_path,
                cfg_deepsomatic.output_path,
                cfg_deepsomatic.model_type,
                output_dir
            )
        }

        if (run_delly2) {
            VARIANT_CALLING_DELLY2_LONGREAD_SOMATIC(
                input_bam_files_ch,
                output_dir,
                uncompressed_fasta,
                cfg_delly2.exclude_tsv_file,
                cfg_delly2.extra_args ?: ''
            )
        }

        if (run_dysgu) {
            VARIANT_CALLING_DYSGU_SOMATIC(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_dysgu.run_extra_args ?: '',
                cfg_dysgu.filter_extra_args ?: '',
                output_dir
            )
        }

        if (run_nanomonsv) {
            VARIANT_CALLING_NANOMONSV(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_nanomonsv.parse_extra_args ?: '',
                cfg_nanomonsv.get_extra_args ?: '',
                output_dir
            )
        }

        if (run_savana) {
            VARIANT_CALLING_SAVANA(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_savana.contigs_txt_file,
                cfg_savana.custom_params_file,
                cfg_savana.run_extra_args ?: '',
                cfg_savana.classify_extra_args ?: '',
                output_dir
            )
        }

        if (run_severus) {
            // Step 1. Phase normal small variant VCF with WhatsHap
            VARIANT_PHASING_WHATSHAP(
                whatshap_input_ch,
                uncompressed_fasta,
                cfg_severus.whatshap_extra_args ?: '',
                output_dir
            )

            // Step 2. Join phased VCF with BAM files for Severus
            // WhatsHap emits: tuple(sample_id, phased_vcf.gz, phased_vcf.gz.tbi)
            // Severus needs: tuple(sample_id, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_vcf)
            severus_input_ch = input_bam_files_ch
                .join(VARIANT_PHASING_WHATSHAP.out.map { sid, vcf, tbi -> tuple(sid, vcf) })
                .map { sid, tbam, tbai, nbam, nbai, vcf -> tuple(sid, tbam, tbai, nbam, nbai, vcf) }

            VARIANT_CALLING_SEVERUS(
                severus_input_ch,
                cfg_severus.vntr_bed_file,
                cfg_severus.extra_args ?: '',
                output_dir
            )
        }

        if (run_svisionpro) {
            VARIANT_CALLING_SVISIONPRO(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_svisionpro.model_file,
                cfg_svisionpro.extra_args ?: '',
                cfg_svisionpro.extract_extra_args ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_LONGREAD_SOMATIC(
        input_bam_files_ch,
        whatshap_input_ch,
        params.reference_genome_fasta_file,
        params.output_dir,
        params.clairs,
        params.deepsomatic,
        params.delly2,
        params.dysgu,
        params.nanomonsv,
        params.savana,
        params.severus,
        params.svisionpro
    )
}

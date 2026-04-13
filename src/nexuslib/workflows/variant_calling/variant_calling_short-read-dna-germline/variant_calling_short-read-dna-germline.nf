#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { decompressFile as decompressFasta }            from '../../../tools/utils'
include { VARIANT_CALLING_CLAIR3 }                       from '../../../subworkflows/variant_calling/variant_calling_clair3/variant_calling_clair3'
include { VARIANT_CALLING_DEEPVARIANT }                  from '../../../subworkflows/variant_calling/variant_calling_deepvariant/variant_calling_deepvariant'
include { VARIANT_CALLING_DELLY2_SHORTREAD_GERMLINE }    from '../../../subworkflows/variant_calling/variant_calling_delly2-sr-germline/variant_calling_delly2-sr-germline'
include { VARIANT_CALLING_DYSGU_GERMLINE }               from '../../../subworkflows/variant_calling/variant_calling_dysgu-germline/variant_calling_dysgu-germline'
include { VARIANT_CALLING_GRIDSS2_GERMLINE }             from '../../../subworkflows/variant_calling/variant_calling_gridss2-germline/variant_calling_gridss2-germline'
include { VARIANT_CALLING_HAPLOTYPECALLER }              from '../../../subworkflows/variant_calling/variant_calling_haplotypecaller/variant_calling_haplotypecaller'
include { VARIANT_CALLING_LUMPY_GERMLINE }               from '../../../subworkflows/variant_calling/variant_calling_lumpy-germline/variant_calling_lumpy-germline'
include { VARIANT_CALLING_MANTA_GERMLINE }               from '../../../subworkflows/variant_calling/variant_calling_manta-germline/variant_calling_manta-germline'
include { VARIANT_CALLING_OCTOPUS_GERMLINE }             from '../../../subworkflows/variant_calling/variant_calling_octopus-germline/variant_calling_octopus-germline'
include { VARIANT_CALLING_PINDEL }                       from '../../../subworkflows/variant_calling/variant_calling_pindel/variant_calling_pindel'
include { VARIANT_CALLING_STRELKA2_GERMLINE }            from '../../../subworkflows/variant_calling/variant_calling_strelka2-germline/variant_calling_strelka2-germline'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         =====================================================================
         Identify germline DNA variants in short-read DNA sequencing BAM files
         =====================================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow variant_calling_short-read-dna-germline.nf -params-file params.yaml [--help]

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
def run_clair3          = run_all || active_methods.contains('clair3')
def run_deepvariant     = run_all || active_methods.contains('deepvariant')
def run_delly2          = run_all || active_methods.contains('delly2')
def run_dysgu           = run_all || active_methods.contains('dysgu')
def run_gridss2         = run_all || active_methods.contains('gridss2')
def run_haplotypecaller = run_all || active_methods.contains('haplotypecaller')
def run_lumpy           = run_all || active_methods.contains('lumpy')
def run_manta           = run_all || active_methods.contains('manta')
def run_octopus         = run_all || active_methods.contains('octopus')
def run_pindel          = run_all || active_methods.contains('pindel')
def run_strelka2        = run_all || active_methods.contains('strelka2')

if (!params.samples_tsv_file)            error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                  error "ERROR: output_dir is required."
if (!params.reference_genome_fasta_file) error "ERROR: reference_genome_fasta_file is required."

if (run_deepvariant) {
    if (!params.deepvariant.input_path)  error "ERROR: deepvariant.input_path is required when running deepvariant."
    if (!params.deepvariant.output_path) error "ERROR: deepvariant.output_path is required when running deepvariant."
}

if (run_delly2) {
    if (!params.delly2.exclude_tsv_file) error "ERROR: delly2.exclude_tsv_file is required when running delly2."
}

if (run_octopus) {
    if (!params.octopus.regions_txt_file) error "ERROR: octopus.regions_txt_file is required when running octopus."
}

def known_methods = [
    'all',
    'clair3',
    'deepvariant',
    'delly2',
    'dysgu',
    'gridss2',
    'haplotypecaller',
    'lumpy',
    'manta',
    'octopus',
    'pindel',
    'strelka2'
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
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_SHORTREAD_GERMLINE {
    take:
        input_bam_files_ch          // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        output_dir
        cfg_clair3
        cfg_deepvariant
        cfg_delly2
        cfg_dysgu
        cfg_gridss2
        cfg_haplotypecaller
        cfg_manta
        cfg_octopus
        cfg_pindel
        cfg_strelka2

    main:
        decompressFasta(reference_genome_fasta_file)
        uncompressed_fasta = decompressFasta.out.f

        if (run_clair3) {
            VARIANT_CALLING_CLAIR3(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_clair3.extra_args ?: '',
                output_dir
            )
        }

        if (run_deepvariant) {
            VARIANT_CALLING_DEEPVARIANT(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_deepvariant.containerization,
                cfg_deepvariant.bin_version,
                cfg_deepvariant.bin_path,
                cfg_deepvariant.input_path,
                cfg_deepvariant.output_path,
                cfg_deepvariant.model_type,
                output_dir
            )
        }

        if (run_delly2) {
            VARIANT_CALLING_DELLY2_SHORTREAD_GERMLINE(
                input_bam_files_ch,
                output_dir,
                uncompressed_fasta,
                cfg_delly2.exclude_tsv_file,
                cfg_delly2.extra_args ?: ''
            )
        }

        if (run_dysgu) {
            VARIANT_CALLING_DYSGU_GERMLINE(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_dysgu.run_extra_args ?: '',
                cfg_dysgu.filter_extra_args ?: '',
                output_dir
            )
        }

        if (run_gridss2) {
            VARIANT_CALLING_GRIDSS2_GERMLINE(
                input_bam_files_ch,
                output_dir,
                uncompressed_fasta,
                cfg_gridss2.extra_args ?: ''
            )
        }

        if (run_haplotypecaller) {
            VARIANT_CALLING_HAPLOTYPECALLER(
                input_bam_files_ch,
                output_dir,
                uncompressed_fasta,
                cfg_haplotypecaller.extra_args ?: '',
                cfg_haplotypecaller.chromosomes
            )
        }

        if (run_lumpy) {
            VARIANT_CALLING_LUMPY_GERMLINE(
                input_bam_files_ch,
                output_dir
            )
        }

        if (run_manta) {
            VARIANT_CALLING_MANTA_GERMLINE(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_manta.config_extra_args ?: '',
                cfg_manta.run_extra_args ?: '',
                output_dir
            )
        }

        if (run_octopus) {
            VARIANT_CALLING_OCTOPUS_GERMLINE(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_octopus.regions_txt_file,
                cfg_octopus.extra_args ?: '',
                output_dir
            )
        }

        if (run_pindel) {
            VARIANT_CALLING_PINDEL(
                input_bam_files_ch,
                output_dir,
                uncompressed_fasta,
                cfg_pindel.extra_args ?: ''
            )
        }

        if (run_strelka2) {
            VARIANT_CALLING_STRELKA2_GERMLINE(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_strelka2.extra_args ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_SHORTREAD_GERMLINE(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.output_dir,
        params.clair3,
        params.deepvariant,
        params.delly2,
        params.dysgu,
        params.gridss2,
        params.haplotypecaller,
        params.manta,
        params.octopus,
        params.pindel,
        params.strelka2
    )
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { VARIANT_CALLING_CLAIRS }                          from '../../../subworkflows/variant_calling/variant_calling_clairs/variant_calling_clairs'
include { VARIANT_CALLING_DEEPSOMATIC }                     from '../../../subworkflows/variant_calling/variant_calling_deepsomatic/variant_calling_deepsomatic'
include { VARIANT_CALLING_DELLY2_SHORTREAD_SOMATIC }        from '../../../subworkflows/variant_calling/variant_calling_delly2-sr-somatic/variant_calling_delly2-sr-somatic'
include { VARIANT_CALLING_DYSGU_SOMATIC }                   from '../../../subworkflows/variant_calling/variant_calling_dysgu-somatic/variant_calling_dysgu-somatic'
include { VARIANT_CALLING_GRIDSS2_SOMATIC }                 from '../../../subworkflows/variant_calling/variant_calling_gridss2-somatic/variant_calling_gridss2-somatic'
include { VARIANT_CALLING_LUMPY_SOMATIC }                   from '../../../subworkflows/variant_calling/variant_calling_lumpy-somatic/variant_calling_lumpy-somatic'
include { VARIANT_CALLING_MANTA_SOMATIC }                   from '../../../subworkflows/variant_calling/variant_calling_manta-somatic/variant_calling_manta-somatic'
include { VARIANT_CALLING_MUTECT2 }                         from '../../../subworkflows/variant_calling/variant_calling_mutect2/variant_calling_mutect2'
include { VARIANT_CALLING_OCTOPUS_SOMATIC }                 from '../../../subworkflows/variant_calling/variant_calling_octopus-somatic/variant_calling_octopus-somatic'
include { VARIANT_CALLING_SEQUENZA }                        from '../../../subworkflows/variant_calling/variant_calling_sequenza/variant_calling_sequenza'
include { VARIANT_CALLING_STRELKA2_SOMATIC }                from '../../../subworkflows/variant_calling/variant_calling_strelka2-somatic/variant_calling_strelka2-somatic'
include { VARIANT_CALLING_SVABA }                           from '../../../subworkflows/variant_calling/variant_calling_svaba/variant_calling_svaba'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ===================================================================
         Identify somatic DNA variants in short-read DNA sequencing BAM files
         ===================================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow variant_calling_short-read-dna-somatic.nf -params-file params.yaml [--help]

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
def run_gridss2         = run_all || active_methods.contains('gridss2')
def run_lumpy           = run_all || active_methods.contains('lumpy')
def run_manta           = run_all || active_methods.contains('manta')
def run_mutect2         = run_all || active_methods.contains('mutect2')
def run_octopus         = run_all || active_methods.contains('octopus')
def run_sequenza        = run_all || active_methods.contains('sequenza')
def run_strelka2        = run_all || active_methods.contains('strelka2')
def run_svaba           = run_all || active_methods.contains('svaba')

if (!params.samples_tsv_file)             error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                   error "ERROR: output_dir is required."
if (!params.reference_genome_fasta_file)  error "ERROR: reference_genome_fasta_file is required."

if (run_deepsomatic) {
    if (!params.deepsomatic.input_path)   error "ERROR: deepsomatic.input_path is required when running deepsomatic."
    if (!params.deepsomatic.output_path)  error "ERROR: deepsomatic.output_path is required when running deepsomatic."
}

if (run_delly2) {
    if (!params.delly2.exclude_tsv_file)   error "ERROR: delly2.exclude_tsv_file is required when running delly2."
}

if (run_octopus) {
    if (!params.octopus.regions_txt_file)    error "ERROR: octopus.regions_txt_file is required when running octopus."
}

if (run_sequenza) {
    if (!params.sequenza.assembly)      error "ERROR: sequenza.assembly is required when running sequenza."
}

def known_methods = [
    'all',
    'clairs',
    'deepsomatic',
    'delly2',
    'dysgu',
    'gridss2',
    'lumpy',
    'manta',
    'mutect2',
    'octopus',
    'sequenza',
    'strelka2',
    'svaba'
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
        "${row.tumor_bam_file}",
        "${row.tumor_bam_bai_file}",
        "${row.normal_bam_file}",
        "${row.normal_bam_bai_file}",
        "${row.tumor_bam_file_no_realign}",
        "${row.tumor_bam_bai_file_no_realign}",
        "${row.normal_bam_file_no_realign}",
        "${row.normal_bam_bai_file_no_realign}",
        "${row.normal_sample_id}") }
    .set { input_bam_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_SHORTREAD_SOMATIC {
    take:
        input_bam_files_ch
        reference_genome_fasta_file
        output_dir
        cfg_clairs
        cfg_deepsomatic
        cfg_delly2
        cfg_dysgu
        cfg_gridss2
        cfg_manta
        cfg_mutect2
        cfg_octopus
        cfg_sequenza
        cfg_strelka2
        cfg_svaba

    main:
        abra2_bam_files_ch_1 = input_bam_files_ch.map { it ->
            tuple(it[0], it[1], it[2], it[3], it[4])
        }
        abra2_bam_files_ch_2 = input_bam_files_ch.map { it ->
            tuple(it[0], it[1], it[2], it[3], it[4], it[9])
        }
        no_realign_bam_files_ch = input_bam_files_ch.map { it ->
            tuple(it[0], it[5], it[6], it[7], it[8])
        }

        if (run_clairs) {
            VARIANT_CALLING_CLAIRS(
                abra2_bam_files_ch_1,
                reference_genome_fasta_file,
                cfg_clairs.extra_args ?: '',
                output_dir
            )
        }

        if (run_deepsomatic) {
            VARIANT_CALLING_DEEPSOMATIC(
                abra2_bam_files_ch_1,
                reference_genome_fasta_file,
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
            VARIANT_CALLING_DELLY2_SHORTREAD_SOMATIC(
                abra2_bam_files_ch_1,
                reference_genome_fasta_file,
                cfg_delly2.exclude_tsv_file,
                cfg_delly2.call_extra_args,
                output_dir
            )
        }

        if (run_dysgu) {
            VARIANT_CALLING_DYSGU_SOMATIC(
                abra2_bam_files_ch_1,
                reference_genome_fasta_file,
                cfg_dysgu.run_extra_args,
                cfg_dysgu.filter_extra_args,
                output_dir
            )
        }

        if (run_gridss2) {
            VARIANT_CALLING_GRIDSS2_SOMATIC(
                no_realign_bam_files_ch,
                reference_genome_fasta_file,
                cfg_gridss2.extra_args,
                cfg_gridss2.somatic_filter_extra_args,
                output_dir
            )
        }

        if (run_lumpy) {
            VARIANT_CALLING_LUMPY_SOMATIC(
                abra2_bam_files_ch_1,
                output_dir
            )
        }

        if (run_manta) {
            VARIANT_CALLING_MANTA_SOMATIC(
                no_realign_bam_files_ch,
                reference_genome_fasta_file,
                cfg_manta.config_extra_args,
                cfg_manta.run_extra_args,
                output_dir
            )
        }

        if (run_mutect2) {
            VARIANT_CALLING_MUTECT2(
                abra2_bam_files_ch_2,
                reference_genome_fasta_file,
                cfg_mutect2.germline_resource_vcf_file,
                cfg_mutect2.panel_of_normals_vcf_file,
                cfg_mutect2.getpileupsummaries_variant_vcf_file,
                cfg_mutect2.extra_args,
                cfg_mutect2.getpileupsummaries_extra_args,
                cfg_mutect2.chromosomes,
                output_dir
            )
        }

        if (run_octopus) {
            VARIANT_CALLING_OCTOPUS_SOMATIC(
                abra2_bam_files_ch_1,
                reference_genome_fasta_file,
                cfg_octopus.regions_txt_file,
                cfg_octopus.extra_args,
                output_dir
            )
        }

        if (run_sequenza) {
            VARIANT_CALLING_SEQUENZA(
                abra2_bam_files_ch_1,
                reference_genome_fasta_file,
                cfg_sequenza.assembly,
                cfg_sequenza.chromosomes,
                cfg_sequenza.sequenzautils_gcwiggle,
                cfg_sequenza.sequenzautils_bam2seqz,
                cfg_sequenza.sequenzautils_seqzbinning,
                output_dir
            )
        }

        if (run_strelka2) {
            VARIANT_CALLING_STRELKA2_SOMATIC(
                abra2_bam_files_ch_1,
                reference_genome_fasta_file,
                cfg_strelka2.extra_args,
                output_dir
            )
        }

        if (run_svaba) {
            VARIANT_CALLING_SVABA(
                abra2_bam_files_ch_1,
                reference_genome_fasta_file,
                cfg_svaba.extra_args,
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_SHORTREAD_SOMATIC(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.output_dir,
        params.clairs,
        params.deepsomatic,
        params.delly2,
        params.dysgu,
        params.gridss2,
        params.manta,
        params.mutect2,
        params.octopus,
        params.sequenza,
        params.strelka2,
        params.svaba
    )
}
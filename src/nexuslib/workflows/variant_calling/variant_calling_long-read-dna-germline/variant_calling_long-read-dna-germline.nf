#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { VARIANT_CALLING_CLAIR3 }                      from '../../../subworkflows/variant_calling/variant_calling_clair3/variant_calling_clair3'
include { VARIANT_CALLING_CUTESV }                      from '../../../subworkflows/variant_calling/variant_calling_cutesv/variant_calling_cutesv'
include { VARIANT_CALLING_DEEPVARIANT }                 from '../../../subworkflows/variant_calling/variant_calling_deepvariant/variant_calling_deepvariant'
include { VARIANT_CALLING_DELLY2_LONGREAD_GERMLINE }    from '../../../subworkflows/variant_calling/variant_calling_delly2-lr-germline/variant_calling_delly2-lr-germline'
include { VARIANT_CALLING_DYSGU_GERMLINE }              from '../../../subworkflows/variant_calling/variant_calling_dysgu-germline/variant_calling_dysgu-germline'
include { VARIANT_CALLING_HIFICNV }                     from '../../../subworkflows/variant_calling/variant_calling_hificnv/variant_calling_hificnv'
include { VARIANT_CALLING_LONGSHOT }                    from '../../../subworkflows/variant_calling/variant_calling_longshot/variant_calling_longshot'
include { VARIANT_CALLING_NANOCALLER }                  from '../../../subworkflows/variant_calling/variant_calling_nanocaller/variant_calling_nanocaller'
include { VARIANT_CALLING_NANOVAR }                     from '../../../subworkflows/variant_calling/variant_calling_nanovar/variant_calling_nanovar'
include { VARIANT_CALLING_PBSV }                        from '../../../subworkflows/variant_calling/variant_calling_pbsv/variant_calling_pbsv'
include { VARIANT_CALLING_SNIFFLES2 }                   from '../../../subworkflows/variant_calling/variant_calling_sniffles2/variant_calling_sniffles2'
include { VARIANT_CALLING_SVIM }                        from '../../../subworkflows/variant_calling/variant_calling_svim/variant_calling_svim'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ====================================================================
         Identify germline DNA variants in long-read DNA sequencing BAM files
         ====================================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow variant_calling_long-read-dna-germline.nf -params-file params.yaml [--help]

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
def run_cutesv          = run_all || active_methods.contains('cutesv')
def run_deepvariant     = run_all || active_methods.contains('deepvariant')
def run_delly2          = run_all || active_methods.contains('delly2')
def run_dysgu           = run_all || active_methods.contains('dysgu')
def run_hificnv         = run_all || active_methods.contains('hificnv')
def run_longshot        = run_all || active_methods.contains('longshot')
def run_nanocaller      = run_all || active_methods.contains('nanocaller')
def run_nanovar         = run_all || active_methods.contains('nanovar')
def run_pbsv            = run_all || active_methods.contains('pbsv')
def run_sniffles2       = run_all || active_methods.contains('sniffles2')
def run_svim            = run_all || active_methods.contains('svim')

if (!params.samples_tsv_file)            error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                  error "ERROR: output_dir is required."
if (!params.reference_genome_fasta_file) error "ERROR: reference_genome_fasta_file is required."

if (run_deepvariant) {
    if (!params.deepvariant.input_path)  error "ERROR: deepvariant.input_path is required when running deepvariant."
    if (!params.deepvariant.output_path) error "ERROR: deepvariant.output_path is required when running deepvariant."
}

if (run_hificnv) {
    if (!params.hificnv.exclude_bed_file) error "ERROR: hificnv.exclude_bed_file is required when running hificnv."
    if (!run_deepvariant) error "ERROR: hificnv requires deepvariant to be enabled (VCF input comes from deepvariant output)."
}

if (run_delly2) {
    if (!params.delly2.exclude_tsv_file) error "ERROR: delly2.exclude_tsv_file is required when running delly2."
}

def known_methods = [
    'all',
    'clair3',
    'cutesv',
    'deepvariant',
    'delly2',
    'dysgu',
    'hificnv',
    'longshot',
    'nanocaller',
    'nanovar',
    'pbsv',
    'sniffles2',
    'svim'
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
workflow VARIANT_CALLING_LONGREAD_GERMLINE {
    take:
        input_bam_files_ch          // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        output_dir
        cfg_clair3
        cfg_cutesv
        cfg_deepvariant
        cfg_delly2
        cfg_dysgu
        cfg_hificnv
        cfg_longshot
        cfg_nanocaller
        cfg_nanovar
        cfg_pbsv
        cfg_sniffles2
        cfg_svim

    main:
        if (run_clair3) {
            VARIANT_CALLING_CLAIR3(
                input_bam_files_ch,
                reference_genome_fasta_file,
                cfg_clair3.extra_args ?: '',
                output_dir
            )
        }

        if (run_cutesv) {
            VARIANT_CALLING_CUTESV(
                input_bam_files_ch,
                reference_genome_fasta_file,
                cfg_cutesv.extra_args ?: '',
                output_dir
            )
        }

        if (run_deepvariant) {
            VARIANT_CALLING_DEEPVARIANT(
                input_bam_files_ch,
                reference_genome_fasta_file,
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
            VARIANT_CALLING_DELLY2_LONGREAD_GERMLINE(
                input_bam_files_ch,
                output_dir,
                reference_genome_fasta_file,
                cfg_delly2.exclude_tsv_file,
                cfg_delly2.extra_args ?: ''
            )
        }

        if (run_dysgu) {
            VARIANT_CALLING_DYSGU_GERMLINE(
                input_bam_files_ch,
                reference_genome_fasta_file,
                cfg_dysgu.run_extra_args ?: '',
                cfg_dysgu.filter_extra_args ?: '',
                output_dir
            )
        }

        if (run_hificnv) {
            // HiFiCNV requires a bgzipped VCF file from DeepVariant
            // DeepVariant emits: tuple(sample_id, vcf.gz, gvcf.gz)
            // Create channel: tuple(sample_id, bam_file, bam_bai_file, vcf_file)
            hificnv_input_ch = input_bam_files_ch
                .join(VARIANT_CALLING_DEEPVARIANT.out.deepvariant_out_ch.map { sid, vcf, gvcf -> tuple(sid, vcf) })

            VARIANT_CALLING_HIFICNV(
                hificnv_input_ch,
                reference_genome_fasta_file,
                cfg_hificnv.exclude_bed_file,
                cfg_hificnv.extra_args ?: '',
                output_dir
            )
        }

        if (run_longshot) {
            VARIANT_CALLING_LONGSHOT(
                input_bam_files_ch,
                reference_genome_fasta_file,
                cfg_longshot.extra_args ?: '',
                output_dir
            )
        }

        if (run_nanocaller) {
            VARIANT_CALLING_NANOCALLER(
                input_bam_files_ch,
                reference_genome_fasta_file,
                cfg_nanocaller.extra_args ?: '',
                output_dir
            )
        }

        if (run_nanovar) {
            VARIANT_CALLING_NANOVAR(
                input_bam_files_ch,
                reference_genome_fasta_file,
                cfg_nanovar.extra_args ?: '',
                output_dir
            )
        }

        if (run_pbsv) {
            VARIANT_CALLING_PBSV(
                input_bam_files_ch,
                reference_genome_fasta_file,
                cfg_pbsv.discover_extra_args ?: '',
                cfg_pbsv.call_extra_args ?: '',
                output_dir
            )
        }

        if (run_sniffles2) {
            VARIANT_CALLING_SNIFFLES2(
                input_bam_files_ch,
                reference_genome_fasta_file,
                cfg_sniffles2.extra_args ?: '',
                output_dir
            )
        }

        if (run_svim) {
            VARIANT_CALLING_SVIM(
                input_bam_files_ch,
                reference_genome_fasta_file,
                cfg_svim.extra_args ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_LONGREAD_GERMLINE(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.output_dir,
        params.clair3,
        params.cutesv,
        params.deepvariant,
        params.delly2,
        params.dysgu,
        params.hificnv,
        params.longshot,
        params.nanocaller,
        params.nanovar,
        params.pbsv,
        params.sniffles2,
        params.svim
    )
}

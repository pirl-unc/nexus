#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { decompressFile as decompressFasta }                    from '../../../tools/utils'
include { VARIANT_CALLING_DEEPVARIANT }                          from '../../../subworkflows/variant_calling/variant_calling_deepvariant/variant_calling_deepvariant'
include { VARIANT_CALLING_LONGSHOT }                             from '../../../subworkflows/variant_calling/variant_calling_longshot/variant_calling_longshot'
include { VARIANT_CALLING_PBSV }                                 from '../../../subworkflows/variant_calling/variant_calling_pbsv/variant_calling_pbsv'

// Aliased phaser imports — DSL2 requires aliasing to call the same subworkflow
// more than once in a single workflow (needed for small_variants_caller='all').
include { HAPLOTYPE_PHASING_HIPHASE as HIPHASE_FROM_DEEPVARIANT }                   from '../../../subworkflows/haplotype_phasing/haplotype_phasing_hiphase/haplotype_phasing_hiphase'
include { HAPLOTYPE_PHASING_HIPHASE as HIPHASE_FROM_LONGSHOT }                      from '../../../subworkflows/haplotype_phasing/haplotype_phasing_hiphase/haplotype_phasing_hiphase'
include { HAPLOTYPE_PHASING_WHATSHAP as WHATSHAP_FROM_DEEPVARIANT }                 from '../../../subworkflows/haplotype_phasing/haplotype_phasing_whatshap/haplotype_phasing_whatshap'
include { HAPLOTYPE_PHASING_WHATSHAP as WHATSHAP_FROM_LONGSHOT }                    from '../../../subworkflows/haplotype_phasing/haplotype_phasing_whatshap/haplotype_phasing_whatshap'
// HapCUT2-WhatsHap intentionally runs ONLY from Longshot's VCF.
// Rationale: HapCUT2's input requirements (strict diploid GTs, no '.' calls,
// alleles in {0,1,2}) make DeepVariant's WGS VCFs incompatible with HapCUT2.
include { HAPLOTYPE_PHASING_HAPCUT2_WHATSHAP as HAPCUT2_WHATSHAP_FROM_LONGSHOT }    from '../../../subworkflows/haplotype_phasing/haplotype_phasing_hapcut2-whatshap/haplotype_phasing_hapcut2-whatshap'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ===============================================================================
         Phase small and structural variants in long-read DNA sequencing BAM files
         (DeepVariant and/or Longshot + pbsv + HiPhase / WhatsHap / HapCUT2-WhatsHap)
         ===============================================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow haplotype_phasing_long-read-dna.nf -params-file params.yaml [--help]

    All parameters are supplied via a params.yaml file. See params.yaml for
    full documentation and defaults.
    """.stripIndent()
    exit 0
}

// ------------------------------------------------------------
// Step 3. Validate inputs
// ------------------------------------------------------------
// Small-variants caller selection. Allowed values:
//   "deepvariant" | "longshot" | "all"
// A selected caller ALWAYS runs and publishes its outputs.
def small_variants_caller = (params.small_variants_caller ?: 'deepvariant').toString().trim().toLowerCase()
def known_callers         = ['deepvariant', 'longshot', 'all']
if (!known_callers.contains(small_variants_caller))
    error "ERROR: small_variants_caller must be one of ${known_callers} (got: '${small_variants_caller}')."

def run_deepvariant     = (small_variants_caller == 'deepvariant' || small_variants_caller == 'all')
// run_longshot_caller — Longshot is selected as a *caller*, so its VCF is
//   chained into the downstream caller-dependent phasers (HiPhase / WhatsHap /
//   HapCUT2-WhatsHap) when those are enabled.
def run_longshot_caller = (small_variants_caller == 'longshot'    || small_variants_caller == 'all')

// Methods (phasers + structural-variant caller).
//   methods: "all"        → enable every phaser method
//   methods: "<list>"     → enable just the listed methods
//   methods: ""  / "none" → no phasers (just whatever caller is selected)
def active_methods       = params.methods.toString().tokenize(',').collect { it.trim().toLowerCase() }
def run_all              = active_methods.contains('all')
def run_pbsv             = run_all || active_methods.contains('pbsv')
def run_hiphase          = run_all || active_methods.contains('hiphase')
def run_whatshap         = run_all || active_methods.contains('whatshap')
def run_hapcut2_whatshap = run_all || active_methods.contains('hapcut2-whatshap')
// run_longshot_phaser — Longshot is selected as a *phaser*, meaning we want
//   its standalone phased VCF + phased BAM published. Crucially, this DOES
//   NOT chain Longshot's VCF into the other phasers — that's controlled by
//   run_longshot_caller (above).
def run_longshot_phaser  = run_all || active_methods.contains('longshot')

// Whether VARIANT_CALLING_LONGSHOT actually runs. True if Longshot is needed
// for either chaining (caller) or standalone publishing (phaser).
def run_longshot = run_longshot_caller || run_longshot_phaser

// Auto-promote dependencies.
if (run_hiphase) {
    if (!run_pbsv) log.info "INFO: HiPhase requested — auto-enabling pbsv (structural variants)."
    run_pbsv = true
}
if (run_longshot_phaser && !run_longshot_caller) {
    log.info "INFO: methods includes 'longshot' — Longshot will run and publish its phased VCF + BAM standalone (its VCF will NOT be chained into other phasers)."
}

if (!params.samples_tsv_file)            error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                  error "ERROR: output_dir is required."
if (!params.reference_genome_fasta_file) error "ERROR: reference_genome_fasta_file is required."

if (run_deepvariant) {
    if (!params.deepvariant.input_path)  error "ERROR: deepvariant.input_path is required when running DeepVariant (small_variants_caller='deepvariant' or 'all')."
    if (!params.deepvariant.output_path) error "ERROR: deepvariant.output_path is required when running DeepVariant (small_variants_caller='deepvariant' or 'all')."
}

def known_methods = ['all', 'pbsv', 'hiphase', 'whatshap', 'hapcut2-whatshap', 'longshot']
active_methods.each { m ->
    if (!known_methods.contains(m)) log.warn "WARNING: unknown method '${m}' — will be ignored. (Use 'small_variants_caller' to pick deepvariant/longshot/all.)"
}

log.info """\
    samples_tsv_file             :   ${params.samples_tsv_file}
    reference_genome_fasta_file  :   ${params.reference_genome_fasta_file}
    output_dir                   :   ${params.output_dir}
    small_variants_caller        :   ${small_variants_caller}
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
workflow HAPLOTYPE_PHASING_LONGREAD_DNA {
    take:
        input_bam_files_ch          // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        output_dir
        cfg_deepvariant
        cfg_longshot
        cfg_pbsv
        cfg_whatshap
        cfg_hapcut2_whatshap

    main:
        decompressFasta(reference_genome_fasta_file)
        uncompressed_fasta = decompressFasta.out.f

        // ----------------------------------------------------
        // Small-variants calling (deepvariant and/or longshot)
        // ----------------------------------------------------
        deepvariant_vcf_ch = Channel.empty()
        longshot_vcf_ch    = Channel.empty()

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
            // DeepVariant emits (sid, vcf.gz, gvcf.gz); keep just the .vcf.gz
            deepvariant_vcf_ch = VARIANT_CALLING_DEEPVARIANT.out.map { sid, vcf, gvcf -> tuple(sid, vcf) }
        }

        if (run_longshot) {
            VARIANT_CALLING_LONGSHOT(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_longshot.extra_args ?: '',
                output_dir
            )
            // Longshot emits (sid, vcf, bam, bai); keep just the vcf
            longshot_vcf_ch = VARIANT_CALLING_LONGSHOT.out.map { sid, vcf, bam, bai -> tuple(sid, vcf) }
        }

        // ----------------------------------------------------
        // Structural-variants calling (pbsv)
        // ----------------------------------------------------
        if (run_pbsv) {
            VARIANT_CALLING_PBSV(
                input_bam_files_ch,
                uncompressed_fasta,
                cfg_pbsv.discover_extra_args ?: '',
                cfg_pbsv.call_extra_args ?: '',
                output_dir
            )
        }

        // ----------------------------------------------------
        // Output layout for caller-dependent phasers.
        //
        // Phaser outputs are nested INSIDE each sample's output folder so that
        // everything for one sample lives under ${output_dir}/${sample_id}/:
        //
        //   ${output_dir}/${sample_id}/
        //     ├── ${sample_id}_deepvariant.vcf.gz                (DeepVariant — caller)
        //     ├── ${sample_id}_longshot.vcf, .bam, .bam.bai      (Longshot    — caller / standalone phaser)
        //     ├── ${sample_id}_pbsv.vcf                          (pbsv        — SV caller)
        //     ├── deepvariant_pbsv_hiphase/                      (HiPhase from DeepVariant)
        //     ├── longshot_pbsv_hiphase/                         (HiPhase from Longshot)
        //     ├── deepvariant_whatshap/                          (WhatsHap from DeepVariant)
        //     ├── longshot_whatshap/                             (WhatsHap from Longshot)
        //     └── longshot_hapcut2-whatshap/                     (HapCUT2-WhatsHap from Longshot — DeepVariant→HapCUT2 intentionally not wired)
        //
        // The subdir name encodes the full pipeline path (caller → SV caller →
        // phaser) so multiple phaser runs against different caller VCFs don't
        // collide. The phaser sub-workflows accept a `subdir` argument that
        // gets appended to ${output_dir}/${sample_id}/${subdir}/.
        // ----------------------------------------------------
        // Default WhatsHap haplotag args, used by both whatshap and hapcut2-whatshap
        // unless cfg_whatshap.haplotag_extra_args overrides.
        whatshap_haplotag_default = '--ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads 4'

        // ----------------------------------------------------
        // HiPhase — run per available small-variants caller's VCF.
        // ----------------------------------------------------
        if (run_hiphase && run_deepvariant) {
            struct_vcf_ch_dv = VARIANT_CALLING_PBSV.out
            HIPHASE_FROM_DEEPVARIANT(
                input_bam_files_ch,
                deepvariant_vcf_ch,
                struct_vcf_ch_dv,
                reference_genome_fasta_file,
                output_dir,
                'deepvariant_pbsv_hiphase'
            )
        }
        if (run_hiphase && run_longshot_caller) {
            struct_vcf_ch_ls = VARIANT_CALLING_PBSV.out
            HIPHASE_FROM_LONGSHOT(
                input_bam_files_ch,
                longshot_vcf_ch,
                struct_vcf_ch_ls,
                reference_genome_fasta_file,
                output_dir,
                'longshot_pbsv_hiphase'
            )
        }

        // ----------------------------------------------------
        // WhatsHap — run per available small-variants caller's VCF.
        // ----------------------------------------------------
        if (run_whatshap && run_deepvariant) {
            whatshap_input_ch_dv = input_bam_files_ch
                .join(deepvariant_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            WHATSHAP_FROM_DEEPVARIANT(
                whatshap_input_ch_dv,
                reference_genome_fasta_file,
                cfg_whatshap.phase_extra_args   ?: '--mapq 20',
                cfg_whatshap.haplotag_extra_args ?: whatshap_haplotag_default,
                output_dir,
                'deepvariant_whatshap'
            )
        }
        if (run_whatshap && run_longshot_caller) {
            whatshap_input_ch_ls = input_bam_files_ch
                .join(longshot_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            WHATSHAP_FROM_LONGSHOT(
                whatshap_input_ch_ls,
                reference_genome_fasta_file,
                cfg_whatshap.phase_extra_args   ?: '--mapq 20',
                cfg_whatshap.haplotag_extra_args ?: whatshap_haplotag_default,
                output_dir,
                'longshot_whatshap'
            )
        }

        // ----------------------------------------------------
        // HapCUT2 + WhatsHap haplotag — runs ONLY from Longshot's VCF.
        // (The DeepVariant→HapCUT2 path is intentionally not wired; see import
        //  comment at the top of this file.)
        // ----------------------------------------------------
        if (run_hapcut2_whatshap && run_longshot_caller) {
            hapcut2_input_ch_ls = input_bam_files_ch
                .join(longshot_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            HAPCUT2_WHATSHAP_FROM_LONGSHOT(
                hapcut2_input_ch_ls,
                reference_genome_fasta_file,
                cfg_hapcut2_whatshap.read_technology              ?: 'pacbio',
                cfg_hapcut2_whatshap.extracthairs_extra_args      ?: '',
                cfg_hapcut2_whatshap.hapcut2_extra_args           ?: '',
                cfg_hapcut2_whatshap.whatshap_haplotag_extra_args ?: whatshap_haplotag_default,
                output_dir,
                'longshot_hapcut2-whatshap'
            )
        } else if (run_hapcut2_whatshap && !run_longshot_caller) {
            log.warn "WARNING: methods includes 'hapcut2-whatshap' but small_variants_caller does not include 'longshot' — HapCUT2 will be skipped (DeepVariant→HapCUT2 path is not wired)."
        }
        // ----------------------------------------------------
        // Longshot phaser — already covered by VARIANT_CALLING_LONGSHOT above.
        // Whenever run_longshot_phaser or run_longshot_caller is true,
        // Longshot runs and publishes its phased VCF + haplotagged BAM
        // standalone. No additional subworkflow invocation is needed for the
        // phaser path; chaining (longshot's VCF -> HiPhase/WhatsHap/HapCUT2-
        // WhatsHap) only happens when run_longshot_caller is true.
        // ----------------------------------------------------
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    HAPLOTYPE_PHASING_LONGREAD_DNA(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.output_dir,
        params.deepvariant,
        params.longshot,
        params.pbsv,
        params.whatshap,
        params.hapcut2_whatshap
    )
}

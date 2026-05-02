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
include { VARIANT_CALLING_CLAIR3 }                               from '../../../subworkflows/variant_calling/variant_calling_clair3/variant_calling_clair3'
include { VARIANT_CALLING_PBSV }                                 from '../../../subworkflows/variant_calling/variant_calling_pbsv/variant_calling_pbsv'

// Aliased phaser imports — DSL2 requires aliasing to call the same subworkflow
// more than once in a single workflow (needed when more than one small-variant
// caller is selected and a phaser must run once per caller's VCF).
// HiPhase runs from DeepVariant or Longshot only — Clair3's VCF column is hard-
// coded to "SAMPLE" rather than the BAM SM tag, which conflicts with pbsv's
// SM-derived column and would fail HiPhase's strict --sample-name check.
include { HAPLOTAGGING_HIPHASE as HIPHASE_FROM_DEEPVARIANT }                   from '../../../subworkflows/haplotagging/haplotagging_hiphase/haplotagging_hiphase'
include { HAPLOTAGGING_HIPHASE as HIPHASE_FROM_LONGSHOT }                      from '../../../subworkflows/haplotagging/haplotagging_hiphase/haplotagging_hiphase'
include { HAPLOTAGGING_WHATSHAP as WHATSHAP_FROM_DEEPVARIANT }                 from '../../../subworkflows/haplotagging/haplotagging_whatshap/haplotagging_whatshap'
include { HAPLOTAGGING_WHATSHAP as WHATSHAP_FROM_LONGSHOT }                    from '../../../subworkflows/haplotagging/haplotagging_whatshap/haplotagging_whatshap'
include { HAPLOTAGGING_WHATSHAP as WHATSHAP_FROM_CLAIR3 }                      from '../../../subworkflows/haplotagging/haplotagging_whatshap/haplotagging_whatshap'
// LongPhase runs from any small-variant caller's VCF (no --sample-name
// constraint like HiPhase, so Clair3→LongPhase works fine).
include { HAPLOTAGGING_LONGPHASE as LONGPHASE_FROM_DEEPVARIANT }               from '../../../subworkflows/haplotagging/haplotagging_longphase/haplotagging_longphase'
include { HAPLOTAGGING_LONGPHASE as LONGPHASE_FROM_LONGSHOT }                  from '../../../subworkflows/haplotagging/haplotagging_longphase/haplotagging_longphase'
include { HAPLOTAGGING_LONGPHASE as LONGPHASE_FROM_CLAIR3 }                    from '../../../subworkflows/haplotagging/haplotagging_longphase/haplotagging_longphase'
// Margin runs from any small-variant caller's VCF. Requires HMM params
// JSON files (one for `margin phase`, one for `margin haplotag`) — the
// user supplies platform-appropriate paths via the margin config block.
include { HAPLOTAGGING_MARGIN as MARGIN_FROM_DEEPVARIANT }                     from '../../../subworkflows/haplotagging/haplotagging_margin/haplotagging_margin'
include { HAPLOTAGGING_MARGIN as MARGIN_FROM_LONGSHOT }                        from '../../../subworkflows/haplotagging/haplotagging_margin/haplotagging_margin'
include { HAPLOTAGGING_MARGIN as MARGIN_FROM_CLAIR3 }                          from '../../../subworkflows/haplotagging/haplotagging_margin/haplotagging_margin'
// HapCUT2-WhatsHap runs from Longshot OR Clair3 — both produce clean diploid GTs.
// Rationale: HapCUT2's input requirements (strict diploid GTs, no '.' calls,
// alleles in {0,1,2}) make DeepVariant's WGS VCFs incompatible with HapCUT2,
// so the DeepVariant→HapCUT2 path is intentionally NOT wired.
include { HAPLOTAGGING_HAPCUT2_WHATSHAP as HAPCUT2_WHATSHAP_FROM_LONGSHOT }    from '../../../subworkflows/haplotagging/haplotagging_hapcut2-whatshap/haplotagging_hapcut2-whatshap'
include { HAPLOTAGGING_HAPCUT2_WHATSHAP as HAPCUT2_WHATSHAP_FROM_CLAIR3 }      from '../../../subworkflows/haplotagging/haplotagging_hapcut2-whatshap/haplotagging_hapcut2-whatshap'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ===============================================================================
         Phase small and structural variants in long-read DNA sequencing BAM files
         (DeepVariant / Longshot / Clair3 + pbsv + HiPhase / WhatsHap / HapCUT2-WhatsHap)
         ===============================================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow haplotagging_long-read-dna.nf -params-file params.yaml [--help]

    All parameters are supplied via a params.yaml file. See params.yaml for
    full documentation and defaults.
    """.stripIndent()
    exit 0
}

// ------------------------------------------------------------
// Step 3. Validate inputs
//
// Three orthogonal selectors (CSV lists, "all", or "none"/""):
//   small_variant_callers       subset of {deepvariant, longshot, clair3}
//   structural_variant_callers  subset of {pbsv}
//   phasing_methods             subset of {whatshap, hiphase, hapcut2-whatshap,
//                                          longphase, margin}
//                               default (when unset/blank): "longphase,margin"
// Cross-axis constraint: hiphase requires pbsv.
// ------------------------------------------------------------

// ---- small_variant_callers ----
def known_svc = ['deepvariant', 'longshot', 'clair3', 'all', 'none']
def svc_raw   = (params.small_variant_callers ?: '').toString().trim().toLowerCase()
def svc_list  = svc_raw.tokenize(',').collect { it.trim() }.findAll { it }
svc_list.each { c ->
    if (!known_svc.contains(c))
        error "ERROR: small_variant_callers contains unknown value '${c}' — allowed: ${known_svc - 'none'} or 'all'."
}
def run_all_svc      = svc_list.contains('all')
def run_deepvariant  = run_all_svc || svc_list.contains('deepvariant')
def run_longshot     = run_all_svc || svc_list.contains('longshot')
def run_clair3       = run_all_svc || svc_list.contains('clair3')
if (!(run_deepvariant || run_longshot || run_clair3))
    error "ERROR: small_variant_callers must include at least one of [deepvariant, longshot, clair3] or 'all' (got: '${svc_raw}')."

// ---- structural_variant_callers ----
def known_svr = ['pbsv', 'all', 'none']
def svr_raw   = (params.structural_variant_callers ?: 'none').toString().trim().toLowerCase()
def svr_list  = svr_raw.tokenize(',').collect { it.trim() }.findAll { it }
svr_list.each { c ->
    if (!known_svr.contains(c))
        error "ERROR: structural_variant_callers contains unknown value '${c}' — allowed: [pbsv], or 'all', or 'none'."
}
// "none" wins if present alongside other values (defensive).
def svr_none    = svr_list.contains('none') || svr_list.isEmpty()
def run_all_svr = svr_list.contains('all') && !svr_none
def run_pbsv    = !svr_none && (run_all_svr || svr_list.contains('pbsv'))

// ---- phasing_methods ----
// Default: run BOTH longphase and margin if the user does not specify
// phasing_methods (or explicitly leaves it blank). Set phasing_methods="none"
// to disable haplotagging entirely.
def known_pm = ['whatshap', 'hiphase', 'hapcut2-whatshap', 'longphase', 'margin', 'all', 'none']
def pm_raw   = (params.phasing_methods ?: 'longphase,margin').toString().trim().toLowerCase()
def pm_list  = pm_raw.tokenize(',').collect { it.trim() }.findAll { it }
pm_list.each { m ->
    if (!known_pm.contains(m))
        log.warn "WARNING: phasing_methods contains unknown value '${m}' — will be ignored."
}
def pm_none              = pm_list.contains('none') || pm_list.isEmpty()
def run_all_pm           = pm_list.contains('all') && !pm_none
def run_whatshap         = !pm_none && (run_all_pm || pm_list.contains('whatshap'))
def run_hiphase          = !pm_none && (run_all_pm || pm_list.contains('hiphase'))
def run_hapcut2_whatshap = !pm_none && (run_all_pm || pm_list.contains('hapcut2-whatshap'))
def run_longphase        = !pm_none && (run_all_pm || pm_list.contains('longphase'))
def run_margin           = !pm_none && (run_all_pm || pm_list.contains('margin'))

// ---- cross-axis constraint ----
if (run_hiphase && !run_pbsv)
    error "ERROR: phasing_methods includes 'hiphase' but structural_variant_callers does not include 'pbsv'. HiPhase phases SNVs+SVs jointly and requires a pbsv VCF — add 'pbsv' to structural_variant_callers, or drop 'hiphase' from phasing_methods."

// ---- haplotag_output ----
def known_ho = ['bam', 'tsv', 'both']
def haplotag_output = (params.haplotag_output ?: 'bam').toString().trim().toLowerCase()
if (!known_ho.contains(haplotag_output))
    error "ERROR: haplotag_output must be one of ${known_ho} (got: '${params.haplotag_output}')."

// ---- required paths ----
if (!params.samples_tsv_file)            error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                  error "ERROR: output_dir is required."
if (!params.reference_genome_fasta_file) error "ERROR: reference_genome_fasta_file is required."

if (run_deepvariant) {
    if (!params.deepvariant?.input_path)  error "ERROR: deepvariant.input_path is required when small_variant_callers includes 'deepvariant'."
    if (!params.deepvariant?.output_path) error "ERROR: deepvariant.output_path is required when small_variant_callers includes 'deepvariant'."
}

log.info """\
    samples_tsv_file             :   ${params.samples_tsv_file}
    reference_genome_fasta_file  :   ${params.reference_genome_fasta_file}
    output_dir                   :   ${params.output_dir}
    small_variant_callers        :   ${svc_raw}  →  [deepvariant=${run_deepvariant}, longshot=${run_longshot}, clair3=${run_clair3}]
    structural_variant_callers   :   ${svr_raw}  →  [pbsv=${run_pbsv}]
    phasing_methods              :   ${pm_raw}  →  [whatshap=${run_whatshap}, hiphase=${run_hiphase}, hapcut2-whatshap=${run_hapcut2_whatshap}, longphase=${run_longphase}, margin=${run_margin}]
    haplotag_output              :   ${haplotag_output}
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
workflow HAPLOTAGGING_LONGREAD_DNA {
    take:
        input_bam_files_ch          // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        output_dir
        cfg_deepvariant
        cfg_longshot
        cfg_clair3
        cfg_pbsv
        cfg_whatshap
        cfg_hapcut2_whatshap
        cfg_longphase
        cfg_margin

    main:
        decompressFasta(reference_genome_fasta_file)
        uncompressed_fasta = decompressFasta.out.f

        // ----------------------------------------------------
        // Small-variants calling (deepvariant and/or longshot and/or clair3)
        // ----------------------------------------------------
        deepvariant_vcf_ch = Channel.empty()
        longshot_vcf_ch    = Channel.empty()
        clair3_vcf_ch      = Channel.empty()

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

        if (run_clair3) {
            VARIANT_CALLING_CLAIR3(
                input_bam_files_ch,
                reference_genome_fasta_file,
                cfg_clair3?.extra_args ?: '--model_path=/opt/models/hifi_sequel2/ --platform=hifi --min_coverage=3',
                output_dir
            )
            // Clair3 emits (sid, "${sid}_clair3_outputs/"); the merged VCF
            // inside is conventionally named merge_output.vcf.gz.
            clair3_vcf_ch = VARIANT_CALLING_CLAIR3.out.map { sid, dir ->
                tuple(sid, file("${dir}/merge_output.vcf.gz"))
            }
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
        //     ├── ${sample_id}_clair3_outputs/                   (Clair3      — caller; merged VCF inside)
        //     ├── ${sample_id}_pbsv.vcf                          (pbsv        — SV caller)
        //     ├── deepvariant_pbsv_hiphase/                      (HiPhase from DeepVariant)
        //     ├── longshot_pbsv_hiphase/                         (HiPhase from Longshot — Clair3→HiPhase intentionally not wired)
        //     ├── deepvariant_whatshap/                          (WhatsHap from DeepVariant)
        //     ├── longshot_whatshap/                             (WhatsHap from Longshot)
        //     ├── clair3_whatshap/                               (WhatsHap from Clair3)
        //     ├── longshot_hapcut2-whatshap/                     (HapCUT2-WhatsHap from Longshot)
        //     ├── clair3_hapcut2-whatshap/                       (HapCUT2-WhatsHap from Clair3 — DeepVariant→HapCUT2 intentionally not wired)
        //     ├── deepvariant_longphase/                         (LongPhase from DeepVariant)
        //     ├── longshot_longphase/                            (LongPhase from Longshot)
        //     ├── clair3_longphase/                              (LongPhase from Clair3)
        //     ├── deepvariant_margin/                            (Margin from DeepVariant)
        //     ├── longshot_margin/                               (Margin from Longshot)
        //     └── clair3_margin/                                 (Margin from Clair3)
        //
        // The subdir name encodes the full pipeline path (caller → SV caller →
        // phaser) so multiple phaser runs against different caller VCFs don't
        // collide. The phaser sub-workflows accept a `subdir` argument that
        // gets appended to ${output_dir}/${sample_id}/${subdir}/.
        // ----------------------------------------------------
        // Default WhatsHap haplotag args, used by both whatshap and hapcut2-whatshap
        // unless cfg_whatshap.haplotag_extra_args overrides.
        whatshap_haplotag_default = '--ignore-read-groups --skip-missing-contigs --output-threads 4'

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
        if (run_hiphase && run_longshot) {
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
        if (run_hiphase && run_clair3 && !run_deepvariant && !run_longshot) {
            log.warn "WARNING: phasing_methods includes 'hiphase' but the only small-variant caller selected is 'clair3'. HiPhase from Clair3 is not wired (Clair3's VCF column 'SAMPLE' conflicts with pbsv's BAM-SM-derived column). HiPhase will be skipped — use 'whatshap' or 'hapcut2-whatshap' for Clair3."
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
        if (run_whatshap && run_longshot) {
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
        if (run_whatshap && run_clair3) {
            whatshap_input_ch_c3 = input_bam_files_ch
                .join(clair3_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            WHATSHAP_FROM_CLAIR3(
                whatshap_input_ch_c3,
                reference_genome_fasta_file,
                cfg_whatshap.phase_extra_args   ?: '--mapq 20',
                cfg_whatshap.haplotag_extra_args ?: whatshap_haplotag_default,
                output_dir,
                'clair3_whatshap'
            )
        }

        // ----------------------------------------------------
        // HapCUT2 + WhatsHap haplotag — runs from Longshot AND/OR Clair3.
        // (The DeepVariant→HapCUT2 path is intentionally not wired; see import
        //  comment at the top of this file.)
        // ----------------------------------------------------
        if (run_hapcut2_whatshap && run_longshot) {
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
        }
        if (run_hapcut2_whatshap && run_clair3) {
            hapcut2_input_ch_c3 = input_bam_files_ch
                .join(clair3_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            HAPCUT2_WHATSHAP_FROM_CLAIR3(
                hapcut2_input_ch_c3,
                reference_genome_fasta_file,
                cfg_hapcut2_whatshap.read_technology              ?: 'pacbio',
                cfg_hapcut2_whatshap.extracthairs_extra_args      ?: '',
                cfg_hapcut2_whatshap.hapcut2_extra_args           ?: '',
                cfg_hapcut2_whatshap.whatshap_haplotag_extra_args ?: whatshap_haplotag_default,
                output_dir,
                'clair3_hapcut2-whatshap'
            )
        }
        if (run_hapcut2_whatshap && !run_longshot && !run_clair3) {
            log.warn "WARNING: phasing_methods includes 'hapcut2-whatshap' but small_variant_callers contains neither 'longshot' nor 'clair3' — HapCUT2 will be skipped (DeepVariant→HapCUT2 path is not wired)."
        }

        // ----------------------------------------------------
        // LongPhase — run per available small-variants caller's VCF.
        // Unlike HiPhase, LongPhase has no --sample-name strict check, so all
        // three caller→LongPhase paths are wired (incl. Clair3→LongPhase).
        //
        // `longphase phase` v2.0.1 requires exactly one of --ont or --pb. Default
        // to --ont when the user did not set phase_extra_args (matches LongPhase's
        // historical default). `longphase haplotag` has no platform flag, so its
        // extra args default to '' (empty).
        // ----------------------------------------------------
        if (run_longphase) {
            // Resolve once so empty-string ('') from the user is preserved
            // (Groovy `?:` would otherwise treat '' as false and substitute
            // the default).
            def lp_phase_args    = (cfg_longphase.containsKey('phase_extra_args')    ? cfg_longphase.phase_extra_args    : '--ont')
            def lp_haplotag_args = (cfg_longphase.containsKey('haplotag_extra_args') ? cfg_longphase.haplotag_extra_args : '')
            if (run_deepvariant) {
                longphase_input_ch_dv = input_bam_files_ch
                    .join(deepvariant_vcf_ch, by: 0)
                    .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
                LONGPHASE_FROM_DEEPVARIANT(
                    longphase_input_ch_dv,
                    reference_genome_fasta_file,
                    lp_phase_args,
                    lp_haplotag_args,
                    output_dir,
                    'deepvariant_longphase'
                )
            }
            if (run_longshot) {
                longphase_input_ch_ls = input_bam_files_ch
                    .join(longshot_vcf_ch, by: 0)
                    .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
                LONGPHASE_FROM_LONGSHOT(
                    longphase_input_ch_ls,
                    reference_genome_fasta_file,
                    lp_phase_args,
                    lp_haplotag_args,
                    output_dir,
                    'longshot_longphase'
                )
            }
            if (run_clair3) {
                longphase_input_ch_c3 = input_bam_files_ch
                    .join(clair3_vcf_ch, by: 0)
                    .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
                LONGPHASE_FROM_CLAIR3(
                    longphase_input_ch_c3,
                    reference_genome_fasta_file,
                    lp_phase_args,
                    lp_haplotag_args,
                    output_dir,
                    'clair3_longphase'
                )
            }
        }

        // ----------------------------------------------------
        // Margin — run per available small-variants caller's VCF.
        // Requires an HMM params JSON file (platform-specific, no default).
        // ----------------------------------------------------
        if (run_margin) {
            def margin_phase_json = (cfg_margin instanceof Map) ? cfg_margin.phase_params_json_file : null
            if (!margin_phase_json) {
                error "ERROR: phasing_methods includes 'margin' but margin.phase_params_json_file is not set. Provide a host-accessible Margin phase HMM parameters JSON (e.g. test/data/indices/margin/phase/allParams.phase_vcf.pb-hifi.json), or drop 'margin' from phasing_methods."
            }
        }
        if (run_margin && run_deepvariant) {
            margin_input_ch_dv = input_bam_files_ch
                .join(deepvariant_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            MARGIN_FROM_DEEPVARIANT(
                margin_input_ch_dv,
                reference_genome_fasta_file,
                cfg_margin.phase_params_json_file,
                cfg_margin.phase_extra_args    ?: '',
                output_dir,
                'deepvariant_margin'
            )
        }
        if (run_margin && run_longshot) {
            margin_input_ch_ls = input_bam_files_ch
                .join(longshot_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            MARGIN_FROM_LONGSHOT(
                margin_input_ch_ls,
                reference_genome_fasta_file,
                cfg_margin.phase_params_json_file,
                cfg_margin.phase_extra_args    ?: '',
                output_dir,
                'longshot_margin'
            )
        }
        if (run_margin && run_clair3) {
            margin_input_ch_c3 = input_bam_files_ch
                .join(clair3_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            MARGIN_FROM_CLAIR3(
                margin_input_ch_c3,
                reference_genome_fasta_file,
                cfg_margin.phase_params_json_file,
                cfg_margin.phase_extra_args    ?: '',
                output_dir,
                'clair3_margin'
            )
        }
        // ----------------------------------------------------
        // Note: when 'longshot' is in small_variant_callers,
        // VARIANT_CALLING_LONGSHOT above already publishes its phased VCF +
        // haplotagged BAM. No standalone-phaser path is needed.
        // ----------------------------------------------------
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    HAPLOTAGGING_LONGREAD_DNA(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.output_dir,
        params.deepvariant,
        params.longshot,
        params.clair3,
        params.pbsv,
        params.whatshap,
        params.hapcut2_whatshap,
        params.longphase ?: [:],
        params.margin ?: [:]
    )
}

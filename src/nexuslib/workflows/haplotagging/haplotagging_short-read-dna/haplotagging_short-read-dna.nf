#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 *
 * Phase small variants in short-read (Illumina) DNA BAM files.
 *
 *   small_variant_callers (one or more):
 *     - DeepVariant            (model_type WGS / WES)
 *     - GATK4 HaplotypeCaller  (per-chromosome -> Picard MergeVcfs)
 *     - Clair3                 (Illumina model)
 *     - Strelka2 (germline)
 *
 *   phasing_methods (zero or more):
 *     - WhatsHap   (phase + haplotag)
 *     - HapCUT2    + WhatsHap haplotag
 *
 * Each caller's small-variants VCF is fed independently into each requested
 * phaser, producing per-(caller, phaser) output folders so multiple
 * combinations don't collide.
 *
 * NOTE: HapCUT2 requires strict diploid genotypes. Some short-read VCFs
 * (especially WGS DeepVariant) may need pre-filtering. If a HapCUT2 task
 * fails on 'Non-diploid VCF entry detected', drop 'hapcut2-whatshap' from
 * phasing_methods for that caller, or pre-filter the VCF.
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { VARIANT_CALLING_DEEPVARIANT }                                             from '../../../subworkflows/variant_calling/variant_calling_deepvariant/variant_calling_deepvariant'
include { VARIANT_CALLING_HAPLOTYPECALLER }                                         from '../../../subworkflows/variant_calling/variant_calling_haplotypecaller/variant_calling_haplotypecaller'
include { VARIANT_CALLING_CLAIR3 }                                                  from '../../../subworkflows/variant_calling/variant_calling_clair3/variant_calling_clair3'
include { VARIANT_CALLING_STRELKA2_GERMLINE }                                       from '../../../subworkflows/variant_calling/variant_calling_strelka2-germline/variant_calling_strelka2-germline'

// Aliased phaser imports - DSL2 requires aliasing to call the same subworkflow
// more than once in a single workflow (one invocation per caller's VCF).
include { HAPLOTAGGING_WHATSHAP as WHATSHAP_FROM_DEEPVARIANT }                 from '../../../subworkflows/haplotagging/haplotagging_whatshap/haplotagging_whatshap'
include { HAPLOTAGGING_WHATSHAP as WHATSHAP_FROM_HAPLOTYPECALLER }             from '../../../subworkflows/haplotagging/haplotagging_whatshap/haplotagging_whatshap'
include { HAPLOTAGGING_WHATSHAP as WHATSHAP_FROM_CLAIR3 }                      from '../../../subworkflows/haplotagging/haplotagging_whatshap/haplotagging_whatshap'
include { HAPLOTAGGING_WHATSHAP as WHATSHAP_FROM_STRELKA2 }                    from '../../../subworkflows/haplotagging/haplotagging_whatshap/haplotagging_whatshap'

include { HAPLOTAGGING_HAPCUT2_WHATSHAP as HAPCUT2_WHATSHAP_FROM_DEEPVARIANT }      from '../../../subworkflows/haplotagging/haplotagging_hapcut2-whatshap/haplotagging_hapcut2-whatshap'
include { HAPLOTAGGING_HAPCUT2_WHATSHAP as HAPCUT2_WHATSHAP_FROM_HAPLOTYPECALLER }  from '../../../subworkflows/haplotagging/haplotagging_hapcut2-whatshap/haplotagging_hapcut2-whatshap'
include { HAPLOTAGGING_HAPCUT2_WHATSHAP as HAPCUT2_WHATSHAP_FROM_CLAIR3 }           from '../../../subworkflows/haplotagging/haplotagging_hapcut2-whatshap/haplotagging_hapcut2-whatshap'
include { HAPLOTAGGING_HAPCUT2_WHATSHAP as HAPCUT2_WHATSHAP_FROM_STRELKA2 }         from '../../../subworkflows/haplotagging/haplotagging_hapcut2-whatshap/haplotagging_hapcut2-whatshap'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ===============================================================================
         Phase small variants in short-read DNA sequencing BAM files
         (DeepVariant / HaplotypeCaller / Clair3 / Strelka2-germline + WhatsHap / HapCUT2)
         ===============================================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow haplotagging_short-read-dna.nf -params-file params.yaml [--help]

    All parameters are supplied via a params.yaml file. See params.yaml for
    full documentation and defaults.
    """.stripIndent()
    exit 0
}

// ------------------------------------------------------------
// Step 3. Validate inputs
// ------------------------------------------------------------

// ---- small_variant_callers ----
def known_svc = ['deepvariant', 'haplotypecaller', 'clair3', 'strelka2-germline', 'all', 'none']
def svc_raw   = (params.small_variant_callers ?: '').toString().trim().toLowerCase()
def svc_list  = svc_raw.tokenize(',').collect { it.trim() }.findAll { it }
svc_list.each { c ->
    if (!known_svc.contains(c))
        error "ERROR: small_variant_callers contains unknown value '${c}' - allowed: ${known_svc - 'none'} or 'all'."
}
def run_all_svc          = svc_list.contains('all')
def run_deepvariant      = run_all_svc || svc_list.contains('deepvariant')
def run_haplotypecaller  = run_all_svc || svc_list.contains('haplotypecaller')
def run_clair3           = run_all_svc || svc_list.contains('clair3')
def run_strelka2         = run_all_svc || svc_list.contains('strelka2-germline')
if (!(run_deepvariant || run_haplotypecaller || run_clair3 || run_strelka2))
    error "ERROR: small_variant_callers must include at least one of [deepvariant, haplotypecaller, clair3, strelka2-germline] or 'all' (got: '${svc_raw}')."

// ---- phasing_methods ----
def known_pm = ['whatshap', 'hapcut2-whatshap', 'all', 'none']
def pm_raw   = (params.phasing_methods ?: '').toString().trim().toLowerCase()
def pm_list  = pm_raw.tokenize(',').collect { it.trim() }.findAll { it }
pm_list.each { m ->
    if (!known_pm.contains(m))
        log.warn "WARNING: phasing_methods contains unknown value '${m}' - will be ignored."
}
def run_all_pm           = pm_list.contains('all')
def pm_none              = pm_list.contains('none') || pm_list.isEmpty()
def run_whatshap         = !pm_none && (run_all_pm || pm_list.contains('whatshap'))
def run_hapcut2_whatshap = !pm_none && (run_all_pm || pm_list.contains('hapcut2-whatshap'))

// ---- haplotag_output ----
def known_ho = ['bam', 'tsv', 'both']
def haplotag_output = (params.haplotag_output ?: 'bam').toString().trim().toLowerCase()
if (!known_ho.contains(haplotag_output))
    error "ERROR: haplotag_output must be one of ${known_ho} (got: '${params.haplotag_output}')."

// ---- core required paths ----
if (!params.samples_tsv_file)            error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                  error "ERROR: output_dir is required."
if (!params.reference_genome_fasta_file) error "ERROR: reference_genome_fasta_file is required."

// ---- caller-specific required paths ----
if (run_deepvariant) {
    if (!params.deepvariant?.input_path)  error "ERROR: deepvariant.input_path is required when small_variant_callers includes 'deepvariant'."
    if (!params.deepvariant?.output_path) error "ERROR: deepvariant.output_path is required when small_variant_callers includes 'deepvariant'."
}

if (run_hapcut2_whatshap) {
    log.warn "NOTE: HapCUT2 requires strict diploid genotypes (alleles in {0,1,2}, no '.' calls). Some short-read VCFs (especially WGS DeepVariant) may need pre-filtering. If a HapCUT2 task fails on 'Non-diploid VCF entry detected', drop 'hapcut2-whatshap' from phasing_methods for that caller."
}

log.info """\
    samples_tsv_file             :   ${params.samples_tsv_file}
    reference_genome_fasta_file  :   ${params.reference_genome_fasta_file}
    output_dir                   :   ${params.output_dir}
    small_variant_callers        :   ${svc_raw}  ->  [deepvariant=${run_deepvariant}, haplotypecaller=${run_haplotypecaller}, clair3=${run_clair3}, strelka2-germline=${run_strelka2}]
    phasing_methods              :   ${pm_raw}  ->  [whatshap=${run_whatshap}, hapcut2-whatshap=${run_hapcut2_whatshap}]
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
// Step 5. Sub-workflow
// ------------------------------------------------------------
workflow HAPLOTAGGING_SHORTREAD_DNA {
    take:
        input_bam_files_ch          // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        output_dir
        cfg_deepvariant
        cfg_haplotypecaller
        cfg_clair3
        cfg_strelka2_germline
        cfg_whatshap
        cfg_hapcut2_whatshap

    main:
        // ----------------------------------------------------
        // Small-variant callers - each emits a (sid, vcf) channel.
        // (Each subworkflow handles its own FASTA decompression/indexing.)
        // ----------------------------------------------------
        deepvariant_vcf_ch     = Channel.empty()
        haplotypecaller_vcf_ch = Channel.empty()
        clair3_vcf_ch          = Channel.empty()
        strelka2_vcf_ch        = Channel.empty()

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
            // DeepVariant emits (sid, vcf.gz, gvcf.gz); keep just the .vcf.gz
            deepvariant_vcf_ch = VARIANT_CALLING_DEEPVARIANT.out.map { sid, vcf, gvcf -> tuple(sid, vcf) }
        }

        if (run_haplotypecaller) {
            VARIANT_CALLING_HAPLOTYPECALLER(
                input_bam_files_ch,
                output_dir,
                reference_genome_fasta_file,
                cfg_haplotypecaller?.extra_args ?: '',
                cfg_haplotypecaller?.chromosomes ?: 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM'
            )
            // Picard MergeVcfs emits (sid, "${sid}_gatk4-haplotypecaller.vcf")
            haplotypecaller_vcf_ch = VARIANT_CALLING_HAPLOTYPECALLER.out
        }

        if (run_clair3) {
            VARIANT_CALLING_CLAIR3(
                input_bam_files_ch,
                reference_genome_fasta_file,
                cfg_clair3?.extra_args ?: '--platform=ilmn --include_all_ctgs',
                output_dir
            )
            // Clair3 emits (sid, "${sid}_clair3_outputs/"); the merged VCF
            // inside is conventionally named merge_output.vcf.gz.
            clair3_vcf_ch = VARIANT_CALLING_CLAIR3.out.map { sid, dir ->
                tuple(sid, file("${dir}/merge_output.vcf.gz"))
            }
        }

        if (run_strelka2) {
            VARIANT_CALLING_STRELKA2_GERMLINE(
                input_bam_files_ch,
                reference_genome_fasta_file,
                cfg_strelka2_germline?.extra_args ?: '',
                output_dir
            )
            // Strelka2-germline emits (sid, "${sid}_strelka2.vcf")
            strelka2_vcf_ch = VARIANT_CALLING_STRELKA2_GERMLINE.out
        }

        // ----------------------------------------------------
        // Output layout for caller-dependent phasers.
        //
        // All outputs land under ${output_dir}/${sample_id}/. Phasers go into
        // their own <caller>_<phaser>/ subfolders so combinations don't collide:
        //
        //   ${output_dir}/${sample_id}/
        //     ├── ${sample_id}_deepvariant.vcf.gz
        //     ├── ${sample_id}_gatk4-haplotypecaller.vcf
        //     ├── ${sample_id}_clair3_outputs/
        //     ├── ${sample_id}_strelka2.vcf
        //     ├── deepvariant_whatshap/
        //     ├── deepvariant_hapcut2-whatshap/
        //     ├── haplotypecaller_whatshap/
        //     ├── haplotypecaller_hapcut2-whatshap/
        //     ├── clair3_whatshap/
        //     ├── clair3_hapcut2-whatshap/
        //     ├── strelka2-germline_whatshap/
        //     └── strelka2-germline_hapcut2-whatshap/
        // ----------------------------------------------------
        whatshap_haplotag_default = '--ignore-read-groups --skip-missing-contigs --output-threads 4'

        // ----------------------------------------------------
        // WhatsHap - one invocation per requested caller.
        // ----------------------------------------------------
        if (run_whatshap && run_deepvariant) {
            ws_in_dv = input_bam_files_ch
                .join(deepvariant_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            WHATSHAP_FROM_DEEPVARIANT(
                ws_in_dv,
                reference_genome_fasta_file,
                cfg_whatshap?.phase_extra_args    ?: '--mapq 20',
                cfg_whatshap?.haplotag_extra_args ?: whatshap_haplotag_default,
                output_dir,
                'deepvariant_whatshap'
            )
        }
        if (run_whatshap && run_haplotypecaller) {
            ws_in_hc = input_bam_files_ch
                .join(haplotypecaller_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            WHATSHAP_FROM_HAPLOTYPECALLER(
                ws_in_hc,
                reference_genome_fasta_file,
                cfg_whatshap?.phase_extra_args    ?: '--mapq 20',
                cfg_whatshap?.haplotag_extra_args ?: whatshap_haplotag_default,
                output_dir,
                'haplotypecaller_whatshap'
            )
        }
        if (run_whatshap && run_clair3) {
            ws_in_c3 = input_bam_files_ch
                .join(clair3_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            WHATSHAP_FROM_CLAIR3(
                ws_in_c3,
                reference_genome_fasta_file,
                cfg_whatshap?.phase_extra_args    ?: '--mapq 20',
                cfg_whatshap?.haplotag_extra_args ?: whatshap_haplotag_default,
                output_dir,
                'clair3_whatshap'
            )
        }
        if (run_whatshap && run_strelka2) {
            ws_in_s2 = input_bam_files_ch
                .join(strelka2_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            WHATSHAP_FROM_STRELKA2(
                ws_in_s2,
                reference_genome_fasta_file,
                cfg_whatshap?.phase_extra_args    ?: '--mapq 20',
                cfg_whatshap?.haplotag_extra_args ?: whatshap_haplotag_default,
                output_dir,
                'strelka2-germline_whatshap'
            )
        }

        // ----------------------------------------------------
        // HapCUT2 + WhatsHap haplotag - one invocation per requested caller.
        //
        // Read technology defaults to 'illumina' for short-read.
        // ----------------------------------------------------
        def read_tech = cfg_hapcut2_whatshap?.read_technology ?: 'illumina'
        if (run_hapcut2_whatshap && run_deepvariant) {
            hc_in_dv = input_bam_files_ch
                .join(deepvariant_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            HAPCUT2_WHATSHAP_FROM_DEEPVARIANT(
                hc_in_dv,
                reference_genome_fasta_file,
                read_tech,
                cfg_hapcut2_whatshap?.extracthairs_extra_args      ?: '',
                cfg_hapcut2_whatshap?.hapcut2_extra_args           ?: '',
                cfg_hapcut2_whatshap?.whatshap_haplotag_extra_args ?: whatshap_haplotag_default,
                output_dir,
                'deepvariant_hapcut2-whatshap'
            )
        }
        if (run_hapcut2_whatshap && run_haplotypecaller) {
            hc_in_hc = input_bam_files_ch
                .join(haplotypecaller_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            HAPCUT2_WHATSHAP_FROM_HAPLOTYPECALLER(
                hc_in_hc,
                reference_genome_fasta_file,
                read_tech,
                cfg_hapcut2_whatshap?.extracthairs_extra_args      ?: '',
                cfg_hapcut2_whatshap?.hapcut2_extra_args           ?: '',
                cfg_hapcut2_whatshap?.whatshap_haplotag_extra_args ?: whatshap_haplotag_default,
                output_dir,
                'haplotypecaller_hapcut2-whatshap'
            )
        }
        if (run_hapcut2_whatshap && run_clair3) {
            hc_in_c3 = input_bam_files_ch
                .join(clair3_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            HAPCUT2_WHATSHAP_FROM_CLAIR3(
                hc_in_c3,
                reference_genome_fasta_file,
                read_tech,
                cfg_hapcut2_whatshap?.extracthairs_extra_args      ?: '',
                cfg_hapcut2_whatshap?.hapcut2_extra_args           ?: '',
                cfg_hapcut2_whatshap?.whatshap_haplotag_extra_args ?: whatshap_haplotag_default,
                output_dir,
                'clair3_hapcut2-whatshap'
            )
        }
        if (run_hapcut2_whatshap && run_strelka2) {
            hc_in_s2 = input_bam_files_ch
                .join(strelka2_vcf_ch, by: 0)
                .map { sid, bam, bai, vcf -> tuple(sid, bam, bai, vcf) }
            HAPCUT2_WHATSHAP_FROM_STRELKA2(
                hc_in_s2,
                reference_genome_fasta_file,
                read_tech,
                cfg_hapcut2_whatshap?.extracthairs_extra_args      ?: '',
                cfg_hapcut2_whatshap?.hapcut2_extra_args           ?: '',
                cfg_hapcut2_whatshap?.whatshap_haplotag_extra_args ?: whatshap_haplotag_default,
                output_dir,
                'strelka2-germline_hapcut2-whatshap'
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    HAPLOTAGGING_SHORTREAD_DNA(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.output_dir,
        params.deepvariant,
        params.haplotypecaller,
        params.clair3,
        params.strelka2_germline,
        params.whatshap,
        params.hapcut2_whatshap
    )
}

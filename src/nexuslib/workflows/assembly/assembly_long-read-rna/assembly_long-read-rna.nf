#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow workflows
// ------------------------------------------------------------
include { ASSEMBLY_ISOFORM }              from '../../../subworkflows/assembly/assembly_isonform/assembly_isonform'
include { ASSEMBLY_ISONFORM_CLUSTERED }   from '../../../subworkflows/assembly/assembly_isonform_clustered/assembly_isonform_clustered'
include { ASSEMBLY_RATTLE }               from '../../../subworkflows/assembly/assembly_rattle/assembly_rattle'
include { ASSEMBLY_RATTLE_CLUSTERED }     from '../../../subworkflows/assembly/assembly_rattle_clustered/assembly_rattle_clustered'
include { ASSEMBLY_RNABLOOM2 }            from '../../../subworkflows/assembly/assembly_rnabloom2/assembly_rnabloom2'
include { ASSEMBLY_RNABLOOM2_CLUSTERED }  from '../../../subworkflows/assembly/assembly_rnabloom2_clustered/assembly_rnabloom2_clustered'
include { ASSEMBLY_STRINGTIE3 }           from '../../../subworkflows/assembly/assembly_stringtie3/assembly_stringtie3'

// ------------------------------------------------------------
// Step 2. Print banner and help
// ------------------------------------------------------------
log.info """\
         ================================================
         Assemble transcripts from long-read RNA samples
         ================================================
         """.stripIndent()

if (params.help) {
    log.info """\
    usage: nexus run --nf-workflow assembly_long-read-rna.nf -params-file params.yaml [--help]

    All parameters are supplied via a params.yaml file. See params.yaml for
    full documentation and defaults.
    """.stripIndent()
    exit 0
}

// ------------------------------------------------------------
// Step 3. Validate inputs
// ------------------------------------------------------------
def active_methods          = params.methods.toString().tokenize(',').collect { it.trim().toLowerCase() }
def run_all                 = active_methods.isEmpty() || active_methods.contains('all')
def run_isonform            = run_all || active_methods.contains('isonform')
def run_isonform_clustered  = run_all || active_methods.contains('isonform_clustered')  || active_methods.contains('isonform-clustered')
def run_rattle              = run_all || active_methods.contains('rattle')
def run_rattle_clustered    = run_all || active_methods.contains('rattle_clustered')    || active_methods.contains('rattle-clustered')
def run_rnabloom2           = run_all || active_methods.contains('rnabloom2')
def run_rnabloom2_clustered = run_all || active_methods.contains('rnabloom2_clustered') || active_methods.contains('rnabloom2-clustered')
def run_stringtie3          = run_all || active_methods.contains('stringtie3')

if (!params.samples_tsv_file)            error "ERROR: samples_tsv_file is required."
if (!params.output_dir)                  error "ERROR: output_dir is required."

if (run_stringtie3) {
    if (!params.reference_genes_gtf_file) error "ERROR: reference_genes_gtf_file is required when running stringtie3."
}

def known_methods = [
    'all',
    'isonform',
    'isonform_clustered',
    'isonform-clustered',
    'rattle',
    'rattle_clustered',
    'rattle-clustered',
    'rnabloom2',
    'rnabloom2_clustered',
    'rnabloom2-clustered',
    'stringtie3'
]
active_methods.each { m ->
    if (!known_methods.contains(m)) log.warn "WARNING: unknown method '${m}' — will be ignored."
}

def run_any_clustered = run_isonform_clustered || run_rattle_clustered || run_rnabloom2_clustered

if (run_any_clustered) {
    def rows = file(params.samples_tsv_file).splitCsv(header: true, sep: '\t')

    def missing = rows.findAll { (it.tsv_file ?: '').trim() && !(it.cluster_method ?: '').trim() }
    if (missing) {
        error "ERROR: cluster_method is required for every row with a tsv_file when a " +
              "*_clustered method is active. Missing for sample_id(s): " +
              "${missing.collect { it.sample_id }.unique().join(', ')}"
    }

    def dups = rows.findAll { (it.tsv_file ?: '').trim() }
                   .groupBy { [it.sample_id, it.cluster_method] }
                   .findAll { k, v -> v.size() > 1 }
                   .keySet()
    if (dups) {
        error "ERROR: duplicate (sample_id, cluster_method) rows would overwrite each " +
              "other: ${dups.collect { it.join('/') }.join(', ')}"
    }

    ['fastq_file', 'bam_file'].each { col ->
        def inconsistent = rows.groupBy { it.sample_id }
                               .findAll { s, rs -> rs.collect { (it[col] ?: '').trim() }.unique().size() > 1 }
                               .keySet()
        if (inconsistent) {
            error "ERROR: rows for the same sample_id disagree on ${col}: " +
                  "${inconsistent.join(', ')}. The de novo and reference-guided methods run " +
                  "once per sample_id, so every row of a sample must repeat the same ${col}."
        }
    }
}

log.info """\
    samples_tsv_file             :   ${params.samples_tsv_file}
    reference_genes_gtf_file     :   ${params.reference_genes_gtf_file}
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
        "${row.fastq_file}",
        "${row.tsv_file}",
        "${row.cluster_method}",
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow ASSEMBLY_LONGREAD_RNA {
    take:
        input_files_ch              // channel: [val(sample_id), path(fastq_file), path(tsv_file), val(cluster_method), path(bam_file), path(bam_bai_file)]
        reference_genes_gtf_file
        output_dir
        cfg_isonform
        cfg_isonform_clustered
        cfg_rattle
        cfg_rattle_clustered
        cfg_rnabloom2
        cfg_rnabloom2_clustered
        cfg_stringtie3

    main:
        fastq_files_ch = input_files_ch.map { it ->
            tuple(it[0], it[1])
        }.unique()
        clustered_files_ch = input_files_ch.map { it ->
            tuple(it[0], it[1], it[2], it[3])
        }
        bam_files_ch = input_files_ch.map { it ->
            tuple(it[0], it[4], it[5])
        }.unique()

        if (run_isonform) {
            ASSEMBLY_ISOFORM(
                fastq_files_ch,
                cfg_isonform.extra_args ?: '',
                output_dir
            )
        }

        if (run_isonform_clustered) {
            ASSEMBLY_ISONFORM_CLUSTERED(
                clustered_files_ch,
                cfg_isonform_clustered.extra_args ?: '',
                output_dir
            )
        }

        if (run_rattle) {
            ASSEMBLY_RATTLE(
                fastq_files_ch,
                cfg_rattle.params_cluster ?: '',
                cfg_rattle.params_correct ?: '',
                cfg_rattle.params_polish ?: '',
                output_dir
            )
        }

        if (run_rattle_clustered) {
            ASSEMBLY_RATTLE_CLUSTERED(
                clustered_files_ch,
                cfg_rattle_clustered.params_cluster ?: '',
                cfg_rattle_clustered.params_correct ?: '',
                cfg_rattle_clustered.params_polish ?: '',
                output_dir
            )
        }

        if (run_rnabloom2) {
            ASSEMBLY_RNABLOOM2(
                fastq_files_ch,
                cfg_rnabloom2.extra_args ?: '',
                output_dir
            )
        }

        if (run_rnabloom2_clustered) {
            ASSEMBLY_RNABLOOM2_CLUSTERED(
                clustered_files_ch,
                cfg_rnabloom2_clustered.extra_args ?: '',
                output_dir
            )
        }

        if (run_stringtie3) {
            ASSEMBLY_STRINGTIE3(
                bam_files_ch,
                reference_genes_gtf_file,
                cfg_stringtie3.extra_args ?: '',
                output_dir
            )
        }
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ASSEMBLY_LONGREAD_RNA(
        input_files_ch,
        params.reference_genes_gtf_file,
        params.output_dir,
        params.isonform,
        params.isonform_clustered,
        params.rattle,
        params.rattle_clustered,
        params.rnabloom2,
        params.rnabloom2_clustered,
        params.stringtie3
    )
}

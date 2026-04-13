#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsSortByName }                  from '../../../tools/samtools'
include { runLonggf }                              from '../../../tools/longgf'
include { decompressFile as decompressGtf }        from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir                   = ''
params.reference_genes_gtf_file     = ''

// Optional arguments
params.params_longgf        = '100 30 100 1 0 3' // <min-overlap-len> <bin_size> <min-map-len> [pseudogene:0(default)/1/other(no filter)] [Secondary_alignment:0(default)] [min_sup_read:2(default)]

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_LONGGF {
    take:
        input_bam_files_ch              // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genes_gtf_file
        params_longgf
        output_dir

    main:
        decompressGtf(reference_genes_gtf_file)

        runSamtoolsSortByName(input_bam_files_ch)
        runLonggf(
            runSamtoolsSortByName.out.f,
            decompressGtf.out.f,
            params_longgf,
            output_dir
        )

    emit:
        runLonggf.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ========================================================================
             Identify RNA variants in long-read RNA sequencing BAM files using Longgf
             ========================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run Longgf.

        usage: nexus run --nf-workflow variant_calling_longgf.nf [required] [optional] [--help]

        required arguments:
            -c                              :   Nextflow .config file.
            -w                              :   Nextflow work directory path.
            --samples_tsv_file              :   TSV file with the following columns: 'sample_id', 'bam_file', 'bam_bai_file'.
            --output_dir                    :   Directory to which output files will be copied.
            --reference_genes_gtf_file      :   Reference genes annotation GTF file.

        optional arguments:
            --params_longgf                 :   Longgf parameters (default: '"100 30 100 1 0 3"').
                                                <min-overlap-len> <bin_size> <min-map-len> [pseudogene:0(default)/1/other(no filter)] [Secondary_alignment:0(default)] [min_sup_read:2(default)]
                                                Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_longgf = (params.params_longgf == true) ? '' : params.params_longgf

    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        reference_genes_gtf_file        :   ${params.reference_genes_gtf_file}
        params_longgf                   :   ${params_longgf}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    VARIANT_CALLING_LONGGF(
        input_bam_files_ch,
        params.reference_genes_gtf_file,
        params_longgf,
        params.output_dir
    )
}

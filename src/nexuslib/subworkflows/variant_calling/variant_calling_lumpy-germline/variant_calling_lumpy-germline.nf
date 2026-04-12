#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runLumpyExpressGermline }  from '../../../tools/lumpy'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_LUMPY_GERMLINE {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        output_dir

    main:
        runLumpyExpressGermline(
            input_bam_files_ch,
            output_dir
        )

    emit:
        runLumpyExpressGermline.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ==================================================================================
             Identify germline variants in paired-end read DNA sequencing BAM files using Lumpy
             ==================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run lumpyexpress (germline mode).

        usage: nexus run --nf-workflow variant_calling_lumpy-germline.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id',
                                        'bam_file',
                                        'bam_bai_file'
            --output_dir            :   Directory to which output files will be copied.
        """.stripIndent()
        exit 0
    }

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    runLumpyExpressGermline(
        input_bam_files_ch,
        params.output_dir
    )
}

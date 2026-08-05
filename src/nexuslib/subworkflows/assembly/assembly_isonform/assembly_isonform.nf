#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runIsonForm }         from '../../../tools/isonform'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                     = ''

// Required arguments
params.samples_tsv_file         = ''
params.output_dir               = ''

// Optional arguments
params.params_ison_pipeline     = '--mode pacbio --iso_abundance 3'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ASSEMBLY_ISOFORM {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        params_ison_pipeline
        output_dir

    main:
        runIsonForm(
            input_fastq_files_ch,
            params_ison_pipeline,
            output_dir
        )

    emit:
        runIsonForm.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =================================================
             Assemble long-read RNA FASTQ files using isONform
             =================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Assemble transcripts using isONform (isON_pipeline.sh).

        usage: nexus run --nf-workflow assembly_isonform.nf [required] [optional] [--help]

        required arguments:
            -c                              :   Nextflow .config file.
            -w                              :   Nextflow work directory path.
            --samples_tsv_file              :   TSV file with the following columns:
                                                'sample_id', 'fastq_file'.
            --output_dir                    :   Directory to which output files will be copied.

        optional arguments:
            --params_ison_pipeline          :   isON_pipeline.sh parameters (default: '"--mode pacbio --iso_abundance 3"').
                                                Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_ison_pipeline = (params.params_ison_pipeline == true) ? '' : params.params_ison_pipeline

    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        params_ison_pipeline            :   ${params_ison_pipeline}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_fastq_files_ch }

    ASSEMBLY_ISOFORM(
        input_fastq_files_ch,
        params_ison_pipeline,
        params.output_dir
    )
}

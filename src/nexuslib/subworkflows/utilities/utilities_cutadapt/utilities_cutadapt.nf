#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runCutadapt }  from '../../../tools/cutadapt'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_cutadapt      = '--poly-a -m 30'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow UTILITIES_CUTADAPT {
    take:
        input_fastq_files_ch             // channel: [val(sample_id), path(fastq_file_1), path(fastq_file_2)]
        params_cutadapt
        output_dir

    main:
        runCutadapt(
            input_fastq_files_ch,
            params_cutadapt,
            output_dir
        )

    emit:
        runCutadapt.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ============
             Run Cutadapt
             ============
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run cutadapt.

        usage: nexus run --nf-workflow utilities_cutadapt.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file_1', 'fastq_file_2'.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_cutadapt       :   cutadapt parameters (default: '--poly-a').
        """.stripIndent()
        exit 0
    }

    def params_cutadapt = (params.params_cutadapt == true) ? '' : params.params_cutadapt

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_cutadapt         :   ${params_cutadapt}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file_1}",
            "${row.fastq_file_2}") }
        .set { input_fastq_files_ch }

    UTILITIES_CUTADAPT(
        input_fastq_files_ch,
        params_cutadapt,
        params.output_dir
    )
}

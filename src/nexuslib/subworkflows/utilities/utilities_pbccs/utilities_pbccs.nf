#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runPbccs }            from '../../../tools/pbccs'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_pbccs         = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow UTILITIES_PBCCS {
    take:
        input_dir_ch             // channel: [val(sample_id), path(rnabloom2_output_dir)]
        params_pbccs
        output_dir

    main:
        runPbccs(
            input_dir_ch,
            params_pbccs,
            output_dir
        )

    emit:
        runPbccs.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =========
             Run PBCCS
             =========
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run pbccs.

        usage: nexus run --nf-workflow utilities_pbccs.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns: 'sample_id', 'subreads_bam_file''.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_pbccs          :   ccs parameters (default: '""').
                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_pbccs = (params.params_pbccs == true) ? '' : params.params_pbccs

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_pbccs            :   ${params_pbccs}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.subreads_bam_file}") }
        .set { input_dir_ch }

    UTILITIES_PBCCS(
        input_dir_ch,
        params_pbccs,
        params.output_dir
    )
}

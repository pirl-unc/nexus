#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runMhcflurry2Predict }    from '../../../tools/mhcflurry2'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir                   = ''

// Optional arguments
params.params_mhcflurry2_predict    = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ANTIGEN_PREDICTION_MHCFLURRY2 {
    take:
        input_csv_files_ch             // channel: [val(sample_id), path(csv_file)]
        params_mhcflurry2_predict
        output_dir

    main:
        runMhcflurry2Predict(
            input_csv_files_ch,
            params_mhcflurry2_predict,
            output_dir
        )

    emit:
        runMhcflurry2Predict.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ================================================
             Predict antigen presentation using MHCflurry 2.0
             ================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run mhcflurry-predict command.

        usage: nexus run --nf-workflow antigen_prediction_mhcflurry2.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'csv_file'. The CSV file should have the columns 'peptide' and 'allele'.
            --output_dir                        :   Directory to which output files will be copied.

        optional arguments:
            --params_mhcflurry2_predict         :   mhcflurry-predict parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_mhcflurry2_predict = (params.params_mhcflurry2_predict == true) ? '' : params.params_mhcflurry2_predict

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        params_mhcflurry2_predict           :   ${params_mhcflurry2_predict}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.csv_file}") }
        .set { input_csv_files_ch }

    ANTIGEN_PREDICTION_MHCFLURRY2(
        input_csv_files_ch,
        params_mhcflurry2_predict,
        params.output_dir
    )
}

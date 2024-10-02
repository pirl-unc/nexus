#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runMhcflurry2 } from '../../modules/mhcflurry2'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.params_mhcflurry2_predict = ''
params.delete_work_dir = false

if (params.params_mhcflurry2_predict == true) {
    params_mhcflurry2_predict = ''
} else {
    params_mhcflurry2_predict = params.params_mhcflurry2_predict
}

// Step 3. Print inputs and help
log.info """\
         ================================================
         Predict antigen presentation using MHCflurry 2.0
         ================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run mhcflurry-predict command.

    usage: nexus run --nf-workflow antigen_presentation_prediction_mhcflurry2.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id',
                                                'mhcflurry2_input_csv_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --params_mhcflurry2_predict         :   mhcflurry-predict parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        params_mhcflurry2_predict           :   ${params_mhcflurry2_predict}
        delete_work_dir                     :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.mhcflurry2_input_csv_file}") }
    .set { input_mhcflurry2_csv_files_ch }

// Step 5. Workflow
workflow ANTIGEN_PRESENTATION_PREDICTION_MHCFLURRY2 {
    take:
        input_mhcflurry2_csv_files_ch             // channel: [val(sample_id), path(mhcflurry2_input_csv_file)]
        params_mhcflurry2_predict
        output_dir

    main:
        runMhcflurry2(
            input_mhcflurry2_csv_files_ch,
            params_mhcflurry2_predict,
            output_dir
        )
}

workflow {
    ANTIGEN_PRESENTATION_PREDICTION_MHCFLURRY2(
        input_mhcflurry2_csv_files_ch,
        params_mhcflurry2_predict,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
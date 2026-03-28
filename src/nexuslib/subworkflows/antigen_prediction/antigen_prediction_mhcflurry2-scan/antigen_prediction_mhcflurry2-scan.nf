#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runMhcflurry2PredictScan }    from '../../../tools/mhcflurry2'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.alleles                          = ''

// Optional arguments
params.params_mhcflurry2_predict_scan   = ''

if (params.params_mhcflurry2_predict_scan == true) {
    params_mhcflurry2_predict_scan = ''
} else {
    params_mhcflurry2_predict_scan = params.params_mhcflurry2_predict_scan
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ================================================
         Predict antigen presentation using MHCflurry 2.0
         ================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run mhcflurry-predict-scan command.

    usage: nexus run --nf-workflow antigen_prediction_mhcflurry2-scan.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'csv_file'.
        --alleles                           :   Alleles to predict.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --params_mhcflurry2_predict_scan    :   mhcflurry-predict-scan parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        alleles                             :   ${params.alleles}
        params_mhcflurry2_predict_scan      :   ${params_mhcflurry2_predict_scan}
    """.stripIndent()
}

// ------------------------------------------------------------
// Step 4. Set channels
// ------------------------------------------------------------
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.csv_file}") }
    .set { input_csv_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow ANTIGEN_PREDICTION_MHCFLURRY2_SCAN {
    take:
        input_csv_files_ch             // channel: [val(sample_id), path(csv_file)]
        alleles
        params_mhcflurry2_predict_scan
        output_dir

    main:
        runMhcflurry2PredictScan(
            input_csv_files_ch,
            alleles,
            params_mhcflurry2_predict_scan,
            output_dir
        )

    emit:
        runMhcflurry2PredictScan.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ANTIGEN_PREDICTION_MHCFLURRY2_SCAN(
        input_csv_files_ch,
        params.alleles,
        params_mhcflurry2_predict_scan,
        params.output_dir
    )
}

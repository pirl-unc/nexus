#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runPbsim3RNA }            from '../../../tools/pbsim3'
include { runSamtoolsMerge }        from '../../../tools/samtools'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir = ''
params.pbsim3_model_file            = ''

// Optional arguments
params.params_pbsim3_mode           = 'errhmm'
params.params_pbsim3                = '--length-min 100 --length-max 1000000 --length-sd 100 --pass-num 10 --strategy trans'

if (params.params_pbsim3_mode == true) {
    params_pbsim3_mode = ''
} else {
    params_pbsim3_mode = params.params_pbsim3_mode
}

if (params.params_pbsim3 == true) {
    params_pbsim3 = ''
} else {
    params_pbsim3 = params.params_pbsim3
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ======================================
         Simulate sequencing reads using PBSIM3
         ======================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run pbsim3.

    usage: nexus run --nf-workflow sequencing_simulation_pbsim3-rna.nf [required] [optional] [--help]

    required arguments:
        -c                          :   Nextflow .config file.
        -w                          :   Nextflow work directory path.
        --samples_tsv_file          :   TSV file with the following columns:
                                        'sample_id', 'transcript_file'.
        --output_dir                :   Directory to which output files will be copied.
        --pbsim3_model_file         :   PBSIM3 error model file (ERRHMM*.model).

    optional arguments:
        --params_pbsim3_mode        :   PBSIM3 mode (default: '"errhmm"'). Either '"errhmm"' or '"qshmm"'.
        --params_pbsim3             :   PBSIM3 parameters (default: "--length-min 100 --length-max 1000000 --pass-num 10 --strategy trans").
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file            :   ${params.samples_tsv_file}
        output_dir                  :   ${params.output_dir}
        pbsim3_model_file           :   ${params.pbsim3_model_file}
        params_pbsim3_mode          :   ${params_pbsim3_mode}
        params_pbsim3               :   ${params_pbsim3}
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
        "${row.transcript_file}") }
    .set { input_transcript_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow SEQUENCING_SIMULATION_PBSIM3_RNA {
    take:
        input_transcript_files_ch            // channel: [val(sample_id), path(transcript_file)]
        pbsim3_model_file
        params_pbsim3_mode
        params_pbsim3
        output_dir

    main:
        // Step 1. Run PBSIM3
        runPbsim3RNA(
            input_transcript_files_ch,
            pbsim3_model_file,
            params_pbsim3_mode,
            params_pbsim3,
            output_dir
        )

    emit:
        runPbsim3RNA.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    SEQUENCING_SIMULATION_PBSIM3_RNA(
        input_transcript_files_ch,
        params.pbsim3_model_file,
        params_pbsim3_mode,
        params_pbsim3,
        params.output_dir
    )
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runRatatoskCorrectFirstPass } from '../../modules/ratatosk'
include { runRatatoskCorrectSecondPass } from '../../modules/ratatosk'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.params_ratatosk_first_correct = ''
params.params_ratatosk_second_correct = ''
// Optional arguments
params.delete_work_dir = false

if (params.params_ratatosk_first_correct == true) {
    params_ratatosk_first_correct = ''
} else {
    params_ratatosk_first_correct = params.params_ratatosk_first_correct
}

if (params.params_ratatosk_second_correct == true) {
    params_ratatosk_second_correct = ''
} else {
    params_ratatosk_second_correct = params.params_ratatosk_second_correct
}

// Step 3. Print inputs and help
log.info """\
         =====================================================================================
         Correct long-read (DNA or RNA) FASTQ files with short-read FASTQ files using Ratatosk
         =====================================================================================

         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Correct reads using Ratatosk (1st pass).
        2. Correct reads using Ratatosk (2nd pass).

    usage: nexus run --nf-workflow long_read_error_correction_ratatosk.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'long_read_fastq_file', 'short_read_r1_fastq_file', 'short_read_r2_fastq_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --params_ratatosk_first_correct     :   Ratatosk first-pass 'correct' parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
        --params_ratatosk_second_correct    :   Ratatosk second-pass 'correct' parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.

    optional arguments:
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        params_ratatosk_first_correct       :   ${params_ratatosk_first_correct}
        params_ratatosk_second_correct      :   ${params_ratatosk_second_correct}
        delete_work_dir                     :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.long_read_fastq_file}",
        "${row.short_read_r1_fastq_file}",
        "${row.short_read_r2_fastq_file}") }
    .set { input_fastq_files_ch }

// Step 5. Workflow
workflow LONG_READ_ERROR_CORRECTION_RATATOSK {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(long_read_fastq_file), path(short_read_r1_fastq_file), path(short_read_r2_fastq_file)]
        params_ratatosk_first_correct
        params_ratatosk_second_correct
        output_dir
    main:
        runRatatoskCorrectFirstPass(
            input_fastq_files_ch,
            params_ratatosk_first_correct
        )
        runRatatoskCorrectSecondPass(
            runRatatoskCorrectFirstPass.out.f,
            params_ratatosk_second_correct,
            output_dir
        )
    emit:
        runRatatoskCorrectSecondPass.out.f
}

workflow {
    LONG_READ_ERROR_CORRECTION_RATATOSK(
        input_fastq_files_ch,
        params_ratatosk_first_correct,
        params_ratatosk_second_correct,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runRatatoskIndexFirstPass } from '../../modules/ratatosk'
include { runRatatoskCorrectFirstPass } from '../../modules/ratatosk'
include { runRatatoskIndexSecondPass } from '../../modules/ratatosk'
include { runRatatoskCorrectSecondPass } from '../../modules/ratatosk'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.ratatosk = ''
params.ratatosk_first_index_params = ''
// Optional arguments
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         =====================================================================================
         Correct long-read (DNA or RNA) FASTQ files with short-read FASTQ files using Ratatosk
         =====================================================================================

         """.stripIndent()

if (params.help) {
    log.info"""\
     workflow:
        1. Index using ratatosk (1st pass).
        2. Correct using ratatosk (1st pass).
        3. Index using ratatosk (2nd pass).
        4. Correct using ratatosk (2nd pass).

    usage: nexus run --nf-workflow long_read_error_correction_ratatosk.nf [required] [optional] [--help]

    required arguments:
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'long_read_fastq_file', 'short_read_r1_fastq_file', 'short_read_r2_fastq_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --ratatosk                          :   ratatosk path.
        --ratatosk_first_index_params       :   ratatosk first-pass 'index' parameters (default: '" "').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
        --ratatosk_first_correct_params     :   ratatosk first-pass 'correct' parameters (default: '" "').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
        --ratatosk_second_index_params      :   ratatosk second-pass 'index' parameters (default: '" "').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
        --ratatosk_second_correct_params    :   ratatosk second-pass 'correct' parameters (default: '" "').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.

    optional arguments:
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        ratatosk                            :   ${params.ratatosk}
        ratatosk_first_index_params         :   ${params.ratatosk_first_index_params}
        ratatosk_first_correct_params       :   ${params.ratatosk_first_correct_params}
        ratatosk_second_index_params        :   ${params.ratatosk_second_index_params}
        ratatosk_second_correct_params      :   ${params.ratatosk_second_correct_params}
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
        ratatosk
        ratatosk_first_index_params
        ratatosk_first_correct_params
        ratatosk_second_index_params
        ratatosk_second_correct_params
        output_dir
    main:
        run_ratatosk_first_index_input_ch = input_fastq_files_ch
        runRatatoskIndexFirstPass(
            run_ratatosk_first_index_input_ch,
            ratatosk,
            ratatosk_first_index_params
        )
        runRatatoskCorrectFirstPass(
            runRatatoskIndexFirstPass.out.f,
            ratatosk,
            ratatosk_first_correct_params
        )
        runRatatoskIndexSecondPass(
            runRatatoskCorrectFirstPass.out.f,
            ratatosk,
            ratatosk_second_index_params
        )
        runRatatoskCorrectSecondPass(
            runRatatoskIndexSecondPass.out.f,
            ratatosk,
            ratatosk_second_correct_params,
            output_dir
        )
    emit:
        runRatatoskCorrectSecondPass.out.f
}

workflow {
    LONG_READ_ERROR_CORRECTION_RATATOSK(
        input_fastq_files_ch,
        params.ratatosk,
        params.ratatosk_first_index_params,
        params.ratatosk_first_correct_params,
        params.ratatosk_second_index_params,
        params.ratatosk_second_correct_params,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

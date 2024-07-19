#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runPicardSingleEndFastqToSam } from '../../modules/picard'
include { runPicardPairedEndFastqToSam } from '../../modules/picard'
include { runSamtoolsSamToBam } from '../../modules/samtools'
include { copyBamFile } from '../../modules/utils'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.read_type = 'single-end'
params.params_picard = ''
params.delete_work_dir = false

if (params.params_picard == true) {
    params_picard = ''
} else {
    params_picard = params.params_picard
}


// Step 3. Print inputs and help
log.info """\
         =================================================
         Convert FASTQ to unaligned BAM files using Picard
         =================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run Picard 'FastqToSam'.
        2. Run samtools to convert SAM files to BAM files.

    usage: nexus run --nf-workflow fastq_to_unaligned_bam.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'fastq_file_1', 'fastq_file_2'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --read_type                         :   read type (default: 'single-end'). Either 'single-end' or 'paired-end'.
        --params_picard                     :   Picard 'FastqToSam' parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        read_type                           :   ${params.read_type}
        params_picard                       :   ${params_picard}
        delete_work_dir                     :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
if (params.read_type == 'single-end') {
    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file_1}") }
        .set { input_fastq_files_ch }
} else {
    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file_1}",
            "${row.fastq_file_2}") }
        .set { input_fastq_files_ch }
}

// Step 5. Workflow
workflow FASTQ_TO_UNALIGNED_BAM {
    take:
        input_fastq_files_ch             // channel: [val(sample_id), path(fastq_file_1), path(fastq_file_2)]
        read_type
        params_picard
        output_dir

    main:
        if (read_type == 'single-end') {
            runPicardSingleEndFastqToSam(
                input_fastq_files_ch,
                params_picard,
                output_dir
            )
            runSamtoolsSamToBam(
                runPicardSingleEndFastqToSam.out.f
            )
        } else {
            runPicardPairedEndFastqToSam(
                input_fastq_files_ch,
                params_picard,
                output_dir
            )
            runSamtoolsSamToBam(
                runPicardPairedEndFastqToSam.out.f
            )
        }
        runSamtoolsSamToBam.out.f.set{ copy_bam_file_input_ch }
        copyBamFile(
            copy_bam_file_input_ch,
            output_dir
        )
}

workflow {
    FASTQ_TO_UNALIGNED_BAM(
        input_fastq_files_ch,
        params.read_type,
        params_picard,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runSamtoolsCoverage } from '../../modules/samtools'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.min_mapping_quality = 20
params.min_base_quality = 20
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         ==============================================
         Compute BAM sequencing coverage using samtools
         ==============================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run samtools coverage.

    usage: nexus run --nf-workflow sequencing_coverage.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --min_mapping_quality               :   Minimum mapping quality (default: 20).
        --min_base_quality                  :   Minimum base quality (default: 20).
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        min_mapping_quality                 :   ${params.min_mapping_quality}
        min_base_quality                    :   ${params.min_base_quality}
        delete_work_dir                     :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow SEQUENCING_COVERAGE {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        min_mapping_quality
        min_base_quality
        output_dir

    main:
        runSamtoolsCoverage(
            input_bam_files_ch,
            min_mapping_quality,
            min_base_quality,
            output_dir
        )
}

workflow {
    SEQUENCING_COVERAGE(
        input_bam_files_ch,
        params.min_mapping_quality,
        params.min_base_quality,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

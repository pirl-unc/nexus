#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runRnaBloom2LongRead } from '../../modules/rnabloom2'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.params_rnabloom2 = '--qual 20 --qual-avg 20 --mincov 3 -ntcard -savebf'
params.delete_work_dir = false

if (params.params_rnabloom2 == true) {
    params_rnabloom2 = ''
} else {
    params_rnabloom2 = params.params_rnabloom2
}

// Step 3. Print inputs and help
log.info """\
         =================================================
         Assemble long-read RNA BAM files using RNA-Bloom2
         =================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Assemble transcripts using RNA-Bloom2.

    usage: nexus run --nf-workflow transcriptome_assembly_rnabloom2.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'fastq_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --params_rnabloom2                  :   RNA-Bloom2 parameters (default: '"--qual 20 --qual-avg 20 --mincov 3 -ntcard -savebf"').
                                                Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        params_rnabloom2                    :   ${params_rnabloom2}
        delete_work_dir                     :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.fastq_file}") }
    .set { input_fastq_files_ch }

// Step 5. Workflow
workflow TRANSCRIPTOME_ASSEMBLY_RNABLOOM2 {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        params_rnabloom2
        output_dir
    main:
        runRnaBloom2LongRead(
            input_fastq_files_ch,
            params_rnabloom2,
            output_dir
        )
    emit:
        runRnaBloom2LongRead.out.f
}

workflow {
    TRANSCRIPTOME_ASSEMBLY_RNABLOOM2(
        input_fastq_files_ch,
        params_rnabloom2,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
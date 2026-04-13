#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runFastQCSingleEndRead }  from '../../../tools/fastqc'
include { runFastQCPairedEndRead }  from '../../../tools/fastqc'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.read_type            = 'single-end'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow UTILITIES_FASTQC {
    take:
        input_fastq_files_ch             // channel: [val(sample_id), path(fastq_file)]
        read_type
        output_dir

    main:
        if (read_type == 'single-end') {
            runFastQCSingleEndRead(
                input_fastq_files_ch,
                output_dir
            )
        } else {
            runFastQCPairedEndRead(
                input_fastq_files_ch,
                output_dir
            )
        }

        fastqc_out_ch = Channel.empty()
        if (read_type == 'single-end') {
            fastqc_out_ch = runFastQCSingleEndRead.out.f
        } else {
            fastqc_out_ch = runFastQCPairedEndRead.out.f
        }

    emit:
        fastqc_out_ch
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ==========
             Run FastQC
             ==========
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run FastQC.

        usage: nexus run --nf-workflow utilities_fastqc.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file_1', 'fastq_file_2'.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --read_type             :   read type (default: 'single-end'). Either 'single-end' or 'paired-end'.
        """.stripIndent()
        exit 0
    }

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        read_type               :   ${params.read_type}
    """.stripIndent()

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

    UTILITIES_FASTQC(
        input_fastq_files_ch,
        params.read_type,
        params.output_dir
    )
}

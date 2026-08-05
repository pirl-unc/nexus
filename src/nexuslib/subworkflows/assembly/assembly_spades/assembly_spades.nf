#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runRNASpades }    from '../../../tools/spades'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_rnaspades     = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ASSEMBLY_SPADES {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file_1), path(fastq_file_2)]
        params_rnaspades
        output_dir

    main:
        runRNASpades(
            input_fastq_files_ch,
            params_rnaspades,
            output_dir
        )

    emit:
        runRNASpades.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ================================================
             Assemble paired-end RNA FASTQ files using Spades
             ================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Assemble transcripts using rnaspades.py.

        usage: nexus run --nf-workflow assembly_spades.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file_1', 'fastq_file_2'.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_rnaspades      :   rnaspades.py parameters (default: '""').
                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_rnaspades = (params.params_rnaspades == true) ? '' : params.params_rnaspades

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_rnaspades        :   ${params_rnaspades}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file_1}",
            "${row.fastq_file_2}") }
        .set { input_fastq_files_ch }

    ASSEMBLY_SPADES(
        input_fastq_files_ch,
        params_rnaspades,
        params.output_dir
    )
}

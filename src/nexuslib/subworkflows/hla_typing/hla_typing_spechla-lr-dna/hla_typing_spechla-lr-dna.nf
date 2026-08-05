#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSpecHLALongRead }    from '../../../tools/spechla'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_spechla       = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow HLA_TYPING_SPECHLA_LONG_READ {
    take:
        input_fastq_files_ch
        params_spechla
        output_dir

    main:
        runSpecHLALongRead(
            input_fastq_files_ch,
            params_spechla,
            output_dir
        )
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ================================================================
             Profile HLA alleles from long-read DNA FASTQ files using SpecHLA
             ================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Profile HLA alleles from paired-end read FASTQ files using SpecHLA.

        usage: nexus run --nf-workflow hla_typing_spechla-lr-dna.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_spechla        :   spechla-long-read extra CLI parameters (default: '""').
                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_spechla = (params.params_spechla == true) ? '' : params.params_spechla

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_spechla          :   ${params_spechla}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_fastq_files_ch }

    HLA_TYPING_SPECHLA_LONG_READ(
        input_fastq_files_ch,
        params_spechla,
        params.output_dir
    )
}

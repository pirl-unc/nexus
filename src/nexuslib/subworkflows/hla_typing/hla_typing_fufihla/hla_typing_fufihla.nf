#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runFufihla }      from '../../../tools/fufihla'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_fufihla       = '--hifi'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow HLA_TYPING_FUFIHLA {
    take:
        input_fastq_files_ch
        params_fufihla
        output_dir

    main:
        runFufihla(
            input_fastq_files_ch,
            params_fufihla,
            output_dir
        )

    emit:
        runFufihla.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ====================================================================
             Profile HLA alleles from long-read DNA FASTQ files using FuFiHLA
             ====================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run fufihla.

        usage: nexus run --nf-workflow hla_typing_fufihla.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_fufihla        :   fufihla extra CLI parameters
                                        (default: '"--hifi"').
                                        Use --hifi for PacBio HiFi reads or --ont for Oxford Nanopore reads.
                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_fufihla = (params.params_fufihla == true) ? '' : params.params_fufihla

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_fufihla          :   ${params_fufihla}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_fastq_files_ch }

    HLA_TYPING_FUFIHLA(
        input_fastq_files_ch,
        params_fufihla,
        params.output_dir
    )
}

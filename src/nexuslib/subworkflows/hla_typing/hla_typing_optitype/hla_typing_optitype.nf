#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runOptiType }  from '../../../tools/optitype'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_optitype      = '--dna'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow HLA_TYPING_OPTITYPE {
    take:
        input_fastq_files_ch
        params_optitype
        output_dir

    main:
        runOptiType(
            input_fastq_files_ch,
            params_optitype,
            output_dir
        )

    emit:
        runOptiType.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =====================================================================================
             Profile HLA alleles using paired-end DNA or RNA sequencing FASTQ files using OptiType
             =====================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Profile HLA alleles using paired-end DNA or RNA sequencing FASTQ files using OptiType.

        usage: nexus run --nf-workflow hla_typing_optitype.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file_1', 'fastq_file_2'.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_optitype       :   OptiTypePipeline.py. parameters (default: '"--dna"').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
        """.stripIndent()
        exit 0
    }

    def params_optitype = (params.params_optitype == true) ? '' : params.params_optitype

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_optitype         :   ${params_optitype}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file_1}",
            "${row.fastq_file_2}") }
        .set { input_fastq_files_ch }

    HLA_TYPING_OPTITYPE(
        input_fastq_files_ch,
        params_optitype,
        params.output_dir
    )
}

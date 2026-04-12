//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runTalon }                               from '../../../tools/talon'
include { decompressFile as decompressGtf }        from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir                   = ''
params.reference_genes_gtf_file     = ''

// Optional arguments
params.params_talon_initdb          = ''
params.params_talon                 = '--create_novel_spliced_genes'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ISOFORM_CHARACTERIZATION_TALON {
    take:
        input_bam_files_ch
        reference_genes_gtf_file
        params_talon_initdb
        params_talon
        output_dir

    main:
        decompressGtf(reference_genes_gtf_file)

        runTalon(
            input_bam_files_ch,
            decompressGtf.out.f,
            params_talon_initdb,
            params_talon,
            output_dir
        )

    emit:
        runTalon.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =================================
             Characterize isoforms using Talon
             =================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run Talon.

        usage: nexus run --nf-workflow isoform_characterization_talon.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'bam_file', 'bam_bai_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genes_gtf_file          :   Reference genes GTF file.

        optional arguments:
            --params_talon_initdb               :   talon_initialize_database parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
            --params_talon                      :   talon parameters (default: '"--create_novel_spliced_genes"').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_talon_initdb = (params.params_talon_initdb == true) ? '' : params.params_talon_initdb
    def params_talon        = (params.params_talon == true) ? '' : params.params_talon

    log.info"""\
        output_dir                          :   ${params.output_dir}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_talon_initdb                 :   ${params_talon_initdb}
        params_talon                        :   ${params_talon}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    ISOFORM_CHARACTERIZATION_TALON(
        input_bam_files_ch,
        params.reference_genes_gtf_file,
        params_talon_initdb,
        params_talon,
        params.output_dir
    )
}

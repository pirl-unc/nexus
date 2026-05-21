#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runBeers2 } from '../../../tools/beers2'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir = ''

// Optional arguments

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow SEQUENCING_SIMULATION_BEERS2 {
    take:
        input_config_files_ch            // channel: [val(sample_id), path(config_file)]
        output_dir

    main:
        runBeers2(
            input_config_files_ch,
            output_dir
        )

    emit:
        runBeers2.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ==================================================================
             Simulate sequencing reads using BEERS2 (Illumina paired-end reads)
             ==================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run BEERS2.

        usage: nexus run --nf-workflow sequencing_simulation_beers2.nf [required] [optional] [--help]

        required arguments:
            -c                          :   Nextflow .config file.
            -w                          :   Nextflow work directory path.
            --samples_tsv_file          :   TSV file with the following columns:
                                            'sample_id', 'config_file'.
            --output_dir                :   Directory to which output files will be copied.

        optional arguments:
        """.stripIndent()
        exit 0
    }

    log.info"""\
        samples_tsv_file            :   ${params.samples_tsv_file}
        output_dir                  :   ${params.output_dir}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.config_file}") }
        .set { input_config_files_ch }

    SEQUENCING_SIMULATION_BEERS2(
        input_config_files_ch,
        params.output_dir
    )
}

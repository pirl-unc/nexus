#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runNeat }  from '../../../tools/neat'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir = ''

// Optional arguments
params.params_neat                = '-R 151 -c 30.0 -p 2 --pe 300 30'
params.min_contig_length          = 1000

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow SEQUENCING_SIMULATION_NEAT {
    take:
        input_fasta_files_ch            // channel: [val(sample_id), path(fasta_file)]
        params_neat
        min_contig_length
        output_dir

    main:
        // Step 1. Run NEAT
        runNeat(
            input_fasta_files_ch,
            params_neat,
            min_contig_length,
            output_dir
        )

    emit:
        runNeat.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ====================================
             Simulate sequencing reads using NEAT
             ====================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run neat.

        usage: nexus run --nf-workflow sequencing_simulation_neat.nf [required] [optional] [--help]

        required arguments:
            -c                          :   Nextflow .config file.
            -w                          :   Nextflow work directory path.
            --samples_tsv_file          :   TSV file with the following columns:
                                            'sample_id', 'fasta_file'.
            --output_dir                :   Directory to which output files will be copied.

        optional arguments:
            --params_neat               :   Neat parameters (default: "-R 151 -c 30.0 -p 2 --pe 300 30").
            --min_contig_length         :   Minimum contig length (default: 1000).
        """.stripIndent()
        exit 0
    }

    def params_neat = (params.params_neat == true) ? '' : params.params_neat

    log.info"""\
        samples_tsv_file            :   ${params.samples_tsv_file}
        output_dir                  :   ${params.output_dir}
        params_neat                 :   ${params_neat}
        min_contig_length           :   ${params.min_contig_length}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fasta_file}") }
        .set { input_fasta_files_ch }

    SEQUENCING_SIMULATION_NEAT(
        input_fasta_files_ch,
        params_neat,
        params.min_contig_length,
        params.output_dir
    )
}

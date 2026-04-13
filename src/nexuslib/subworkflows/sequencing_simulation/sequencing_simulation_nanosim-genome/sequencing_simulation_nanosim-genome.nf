#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runNanoSimGenome }        from '../../../tools/nanosim'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir                   = ''
params.model_prefix                 = ''

// Optional arguments
params.params_nanosim               = '--coverage 30 --homopolymer --KmerBias 5'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow SEQUENCING_SIMULATION_NANOSIM {
    take:
        input_fasta_files_ch            // channel: [val(sample_id), path(fasta_file)]
        model_prefix
        params_nanosim
        output_dir

    main:
        runNanoSimGenome(
            input_fasta_files_ch,
            model_prefix,
            params_nanosim,
            output_dir
        )

    emit:
        runNanoSimGenome.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =======================================
             Simulate sequencing reads using NanoSim
             =======================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run pbsim3.

        usage: nexus run --nf-workflow sequencing_simulation_nanosim-genome.nf [required] [optional] [--help]

        required arguments:
            -c                          :   Nextflow .config file.
            -w                          :   Nextflow work directory path.
            --samples_tsv_file          :   TSV file with the following columns:
                                            'sample_id', 'fasta_file'.
            --output_dir                :   Directory to which output files will be copied.
            --model_prefix              :   NanoSim model prefix.

        optional arguments:
            --params_nanosim            :   simulator.py genome parameters (default: "--coverage 30 --homopolymer").
        """.stripIndent()
        exit 0
    }

    def params_nanosim = (params.params_nanosim == true) ? '' : params.params_nanosim

    log.info"""\
        samples_tsv_file            :   ${params.samples_tsv_file}
        output_dir                  :   ${params.output_dir}
        model_prefix                :   ${params.model_prefix}
        params_nanosim              :   ${params_nanosim}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fasta_file}") }
        .set { input_fasta_files_ch }

    SEQUENCING_SIMULATION_NANOSIM(
        input_fasta_files_ch,
        params.model_prefix,
        params_nanosim,
        params.output_dir
    )
}

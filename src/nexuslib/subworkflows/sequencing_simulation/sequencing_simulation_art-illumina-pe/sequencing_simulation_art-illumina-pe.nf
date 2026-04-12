#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runArtIlluminaPairedEnd } from '../../../tools/art'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir = ''

// Optional arguments
params.params_art                   = '-l 150 -f 15 -ss HSXt -m 200 -s 10 -na'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow SEQUENCING_SIMULATION_ART_ILLUMINA_PE {
    take:
        input_fasta_files_ch            // channel: [val(sample_id), path(fasta_file)]
        params_art
        output_dir

    main:
        runArtIlluminaPairedEnd(
            input_fasta_files_ch,
            params_art,
            output_dir
        )

    emit:
        runArtIlluminaPairedEnd.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===============================================================
             Simulate sequencing reads using ART (Illumina paired-end reads)
             ===============================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run ART (art_illumina -p).

        usage: nexus run --nf-workflow sequencing_simulation_art-illumina-pe.nf [required] [optional] [--help]

        required arguments:
            -c                          :   Nextflow .config file.
            -w                          :   Nextflow work directory path.
            --samples_tsv_file          :   TSV file with the following columns:
                                            'sample_id', 'fasta_file'.
            --output_dir                :   Directory to which output files will be copied.

        optional arguments:
            --params_art                :   art_illumina parameters (default: "-l 151 -f 15 -ss HSXt -m 200 -s 10").
        """.stripIndent()
        exit 0
    }

    def params_art = (params.params_art == true) ? '' : params.params_art

    log.info"""\
        samples_tsv_file            :   ${params.samples_tsv_file}
        output_dir                  :   ${params.output_dir}
        params_art                  :   ${params_art}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fasta_file}") }
        .set { input_fasta_files_ch }

    SEQUENCING_SIMULATION_ART_ILLUMINA_PE(
        input_fasta_files_ch,
        params_art,
        params.output_dir
    )
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runGeluster }    from '../../../tools/geluster'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_geluster      = '--rform fq --seqType PacBio'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow READ_CLUSTERING_GELUSTER {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        params_geluster
        output_dir

    main:
        runGeluster(
            input_fastq_files_ch,
            params_geluster,
            output_dir
        )

    emit:
        runGeluster.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ==========================================
             Cluster long-read RNA reads using GeLuster
             ==========================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Cluster sequencing reads using GeLuster.

        usage: nexus run --nf-workflow read_clustering_geluster.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_geluster       :   GeLuster parameters (default: '"--seqType PacBio"').
                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_geluster   = (params.params_geluster == true) ? '' : params.params_geluster

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_geluster       :   ${params_geluster}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_files_ch }

    READ_CLUSTERING_GELUSTER(
        input_files_ch,
        params_geluster,
        params.output_dir
    )
}

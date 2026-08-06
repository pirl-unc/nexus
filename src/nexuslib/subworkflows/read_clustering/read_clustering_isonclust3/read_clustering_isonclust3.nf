#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runIsonClust3 }    from '../../../tools/isonclust3'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_isonclust3    = '--mode pacbio --no-fastq --n 3'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow READ_CLUSTERING_ISONCLUST3 {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        params_isonclust3
        output_dir

    main:
        runIsonClust3(
            input_fastq_files_ch,
            params_isonclust3,
            output_dir
        )

    emit:
        f   = runIsonClust3.out.f
        tsv = runIsonClust3.out.tsv
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ============================================
             Cluster long-read RNA reads using isONclust3
             ============================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Cluster sequencing reads using isONclust3.

        usage: nexus run --nf-workflow read_clustering_isonclust3.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_isonclust3     :   isONclust3 parameters (default: '"--no-fastq --n 3"').
                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_isonclust3   = (params.params_isonclust3 == true) ? '' : params.params_isonclust3

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_isonclust3       :   ${params_isonclust3}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_files_ch }

    READ_CLUSTERING_ISONCLUST3(
        input_files_ch,
        params_isonclust3,
        params.output_dir
    )
}

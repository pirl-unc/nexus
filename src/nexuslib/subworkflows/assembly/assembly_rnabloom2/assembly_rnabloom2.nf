#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runRnaBloom2LongRead }    from '../../../tools/rnabloom2'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_rnabloom2     = '--qual 20 --qual-avg 20 --mincov 3 -ntcard -savebf -chimera -lrpb'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ASSEMBLY_RNABLOOM2 {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        params_rnabloom2
        output_dir

    main:
        runRnaBloom2LongRead(
            input_fastq_files_ch,
            params_rnabloom2,
            output_dir
        )

    emit:
        runRnaBloom2LongRead.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===================================================
             Assemble long-read RNA FASTQ files using RNA-Bloom2
             ===================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Assemble transcripts using RNA-Bloom2.

        usage: nexus run --nf-workflow assembly_rnabloom2.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_rnabloom2      :   RNA-Bloom2 parameters (default: '"--qual 20 --qual-avg 20 --mincov 3 -ntcard -savebf -chimera -lrpb"').
                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_rnabloom2 = (params.params_rnabloom2 == true) ? '' : params.params_rnabloom2

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_rnabloom2        :   ${params_rnabloom2}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_fastq_files_ch }

    ASSEMBLY_RNABLOOM2(
        input_fastq_files_ch,
        params_rnabloom2,
        params.output_dir
    )
}

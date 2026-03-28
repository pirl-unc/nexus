//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runIsoseqCluster }    from '../../../tools/isoseq'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ==================================
         Characterize isoforms using Isoseq
         ==================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run isoseq 'cluster2'.

    usage: nexus run --nf-workflow isoform_characterization_isoseq.nf [required] [optional] [--help]

    required arguments:
        -c                      :   Nextflow .config file.
        -w                      :   Nextflow work directory path.
        --samples_tsv_file      :   TSV file with the following columns:
                                    'sample_id', 'bam_file'.
        --output_dir            :   Directory to which output files will be copied.

    optional arguments:
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
    """.stripIndent()
}

// ------------------------------------------------------------
// Step 4. Set channels
// ------------------------------------------------------------
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.bam_file}") }
    .set { input_bam_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow ISOFORM_CHARACTERIZATION_ISOSEQ {
    take:
        input_bam_files_ch            // channel: [val(sample_id), path(bam_file)]
        output_dir

    main:
        runIsoseqCluster(
            input_bam_files_ch,
            output_dir
        )

    emit:
        runIsoseqCluster.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ISOFORM_CHARACTERIZATION_ISOSEQ(
        input_bam_files_ch,
        params.output_dir
    )
}

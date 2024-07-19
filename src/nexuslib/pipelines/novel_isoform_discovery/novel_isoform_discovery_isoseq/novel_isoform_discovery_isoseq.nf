//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runIsoseqCluster } from '../../modules/isoseq'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.output_dir = ''
// Optional arguments
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         ====================================
         Novel isoform discovery using Isoseq
         ====================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run isoseq 'cluster2'.

    usage: nexus run --nf-workflow novel_isoform_discovery_isoseq.nf [required] [optional] [--help]

    required arguments:
        -c                              :   Nextflow .config file.
        -w                              :   Nextflow work directory path.
        --samples_tsv_file              :   TSV file with the following columns:
                                            'sample_id', 'bam_file'.
        --output_dir                    :   Directory to which output files will be copied.

    optional arguments:
        --delete_work_dir               :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        delete_work_dir                 :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.bam_file}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow NOVEL_ISOFORM_DISCOVERY_ISOSEQ {
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

workflow {
    NOVEL_ISOFORM_DISCOVERY_ISOSEQ(
        input_bam_files_ch,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runLumpyExpressTumorNormal } from '../../modules/lumpy'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         =================================================================================
         Identify somatic variants in paired-end read DNA sequencing BAM files using Lumpy
         =================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run lumpyexpress (tumor and normal mode).

    usage: nexus run --nf-workflow paired-end_read_dna_variant_calling_lumpy.nf [required] [optional] [--help]

    required arguments:
        -c                                          :   Nextflow .config file.
        -w                                          :   Nextflow work directory path.
        --samples_tsv_file                          :   TSV file with the following columns:
                                                        'sample_id',
                                                        'tumor_bam_file',
                                                        'tumor_bam_bai_file',
                                                        'normal_bam_file',
                                                        'normal_bam_bai_file'
        --output_dir                                :   Directory to which output files will be copied.

    optional arguments:
        --delete_work_dir                           :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                            :   ${params.samples_tsv_file}
        output_dir                                  :   ${params.output_dir}
        delete_work_dir                             :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.tumor_bam_file}",
        "${row.tumor_bam_bai_file}",
        "${row.normal_bam_file}",
        "${row.normal_bam_bai_file}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow PAIRED_END_READ_DNA_VARIANT_CALLING_LUMPY {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(tumor_sample_id), val(normal_sample_id)]
        output_dir

    main:
        runLumpyExpressTumorNormal(
            input_bam_files_ch,
            output_dir
        )
}

workflow {
    PAIRED_END_READ_DNA_VARIANT_CALLING_LUMPY(
        input_bam_files_ch,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
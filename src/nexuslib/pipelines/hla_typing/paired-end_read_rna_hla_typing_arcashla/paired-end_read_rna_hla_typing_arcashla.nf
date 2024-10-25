#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runArcasHlaPairedEndMode } from '../../modules/arcashla'

// Step 2. Input parameters
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         ============================================================================
         Profile HLA alleles using paired-end RNA sequencing BAM files using arcasHLA
         ============================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Profile HLA alleles using paired-end RNA sequencing BAM files using arcasHLA.

    usage: nexus run --nf-workflow paired-end_read_rna_hla_typing_arcashla.nf [required] [optional] [--help]

    required arguments:
        -c                              :   Nextflow .config file.
        -w                              :   Nextflow work directory path.
        --samples_tsv_file              :   TSV file with the following columns:
                                            'sample_id', 'bam_file', 'bam_bai_file'.
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
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow PAIRED_END_RNA_READ_HLA_TYPING_ARCASHLA {
    take:
        input_bam_files_ch
        output_dir
    main:
        runArcasHlaPairedEndMode(
            input_bam_files_ch,
            output_dir
        )
    emit:
        runArcasHlaPairedEndMode.out.f
}

workflow {
    PAIRED_END_RNA_READ_HLA_TYPING_ARCASHLA(
        input_bam_files_ch,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
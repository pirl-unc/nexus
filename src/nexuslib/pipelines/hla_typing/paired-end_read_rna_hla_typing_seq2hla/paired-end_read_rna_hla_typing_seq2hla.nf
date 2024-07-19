#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runSeq2HLA } from '../../modules/seq2hla'

// Step 2. Input parameters
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.params_seq2hla = ''
params.delete_work_dir = false

if (params.params_seq2hla == true) {
    params_seq2hla = ''
} else {
    params_seq2hla = params.params_seq2hla
}

// Step 3. Print inputs and help
log.info """\
         =============================================================================
         Profile HLA alleles using paired-end RNA sequencing FASTQ files using Seq2HLA
         =============================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Profile HLA alleles using paired-end RNA sequencing FASTQ files using Seq2HLA.

    usage: nexus run --nf-workflow paired-end_read_rna_hla_typing_seq2hla.nf [required] [optional] [--help]

    required arguments:
        -c                              :   Nextflow .config file.
        -w                              :   Nextflow work directory path.
        --samples_tsv_file              :   TSV file with the following columns:
                                            'sample_id', 'fastq_file_1', 'fastq_file_2'.
        --output_dir                    :   Directory to which output files will be copied.

    optional arguments:
        --params_seq2hla                :   Seq2HLA parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
        --delete_work_dir               :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        params_seq2hla                  :   ${params.params_seq2hla}
        delete_work_dir                 :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.fastq_file_1}",
        "${row.fastq_file_2}") }
    .set { input_fastq_files_ch }

// Step 5. Workflow
workflow PAIRED_END_RNA_READ_HLA_TYPING_SEQ2HLA {
    take:
        input_fastq_files_ch
        params_seq2hla
        output_dir
    main:
        runSeq2HLA(
            input_fastq_files_ch,
            params_seq2hla,
            output_dir
        )
    emit:
        runSeq2HLA.out.f
}

workflow {
    PAIRED_END_RNA_READ_HLA_TYPING_SEQ2HLA(
        input_fastq_files_ch,
        params.params_seq2hla,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
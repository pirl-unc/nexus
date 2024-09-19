#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runKallistoQuantPairedEndReads } from '../../modules/kallisto'

// Step 2. Input parameters
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.kallisto_index_file = ''
// Optional arguments
params.params_kallisto_quant = ''
params.delete_work_dir = false

if (params.params_kallisto_quant == true) {
    params_kallisto_quant = ''
} else {
    params_kallisto_quant = params.params_kallisto_quant
}

// Step 3. Print inputs and help
log.info """\
         ===========================================================
         Quantify RNA in pairend-end read FASTQ files using Kallisto
         ===========================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Quantify RNA in paired-end read FASTQ files using Kallisto.

    usage: nexus run --nf-workflow paired-end_read_rna_quantification_kallisto.nf [required] [optional] [--help]

    required arguments:
        -c                          :   Nextflow .config file.
        -w                          :   Nextflow work directory path.
        --samples_tsv_file          :   TSV file with the following columns:
                                        'sample_id', 'fastq_file_1', 'fastq_file_2'.
        --kallisto_index_file       :   Kallisto index file.
        --output_dir                :   Directory to which output files will be copied.

    optional arguments:
        --params_kallisto_quant     :   Kallisto quant parameters (default: '""').
                                        Note that the parameters need to be wrapped in quotes.
        --delete_work_dir           :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file            :   ${params.samples_tsv_file}
        kallisto_index_file         :   ${params.kallisto_index_file}
        params_kallisto_quant       :   ${params_kallisto_quant}
        output_dir                  :   ${params.output_dir}
        delete_work_dir             :   ${params.delete_work_dir}
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
workflow PAIRED_END_READ_RNA_QUANTIFICATION_KALLISTO {
    take:
        input_fastq_files_ch
        kallisto_index_file
        params_kallisto_quant
        output_dir
    main:
        runKallistoQuantPairedEndReads(
            input_fastq_files_ch,
            kallisto_index_file,
            params_kallisto_quant,
            output_dir
        )
    emit:
        runKallistoQuantPairedEndReads.out.f
}

workflow {
    PAIRED_END_READ_RNA_QUANTIFICATION_KALLISTO(
        input_fastq_files_ch,
        params.kallisto_index_file,
        params_kallisto_quant,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
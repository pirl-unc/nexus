#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runTranSigner } from '../../modules/transigner'

// Step 2. Input parameters
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.reference_transcriptome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.transcripts.fa'
params.params_transigner_align = ''
params.params_transigner_prefilter = '--filter -tp -500 -fp -600'
params.params_transigner_em = '--drop --use-score'
// Optional arguments
params.delete_work_dir = false

if (params.params_transigner_align == true) {
    params_transigner_align = ''
} else {
    params_transigner_align = params.params_transigner_align
}

if (params.params_transigner_prefilter == true) {
    params_transigner_prefilter = ''
} else {
    params_transigner_prefilter = params.params_transigner_prefilter
}

if (params.params_transigner_em == true) {
    params_transigner_em = ''
} else {
    params_transigner_em = params.params_transigner_em
}

// Step 3. Print inputs and help
log.info """\
         ======================================================
         Quantify RNA in long-read FASTQ files using TranSigner
         ======================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Quantify RNA in long-read FASTQ files using TranSigner.

    usage: nexus run --nf-workflow long_read_rna_quantification_transigner.nf [required] [optional] [--help]

    required arguments:
        -c                                      :   Nextflow .config file.
        -w                                      :   Nextflow work directory path.
        --samples_tsv_file                      :   TSV file with the following columns:
                                                    'sample_id', 'fastq_file'.
        --reference_transcriptome_fasta_file    :   Reference transcriptome FASTA file
                                                    (default: '/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.transcripts.fa').
        --params_transigner_align               :   TranSigner align parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes
                                                    and a space at the end of the string is necessary.
        --params_transigner_prefilter           :   TranSigner prefilter parameters (default: '"--filter -tp -500 -fp -600"').
                                                    Note that the parameters need to be wrapped in quotes
                                                    and a space at the end of the string is necessary.
        --params_transigner_em                  :   TranSigner em parameters (default: '"--drop --use-score"').
                                                    Note that the parameters need to be wrapped in quotes
                                                    and a space at the end of the string is necessary.
        --output_dir                            :   Directory to which output files will be copied.

    optional arguments:
        --delete_work_dir                       :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                        :   ${params.samples_tsv_file}
        reference_transcriptome_fasta_file      :   ${params.reference_transcriptome_fasta_file}
        params_transigner_align                 :   ${params_transigner_align}
        params_transigner_prefilter             :   ${params_transigner_prefilter}
        params_transigner_em                    :   ${params_transigner_em}
        output_dir                              :   ${params.output_dir}
        delete_work_dir                         :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.fastq_file}") }
    .set { input_fastq_files_ch }

// Step 5. Workflow
workflow LONG_READ_RNA_QUANTIFICATION_TRANSIGNER {
    take:
        input_fastq_files_ch
        reference_transcriptome_fasta_file
        params_transigner_align
        params_transigner_prefilter
        params_transigner_em
        output_dir
    main:
        runTranSigner(
            input_fastq_files_ch,
            reference_transcriptome_fasta_file,
            params_transigner_align,
            params_transigner_prefilter,
            params_transigner_em,
            output_dir
        )
    emit:
        runTranSigner.out.f
}

workflow {
    LONG_READ_RNA_QUANTIFICATION_TRANSIGNER(
        input_fastq_files_ch,
        params.reference_transcriptome_fasta_file,
        params_transigner_align,
        params_transigner_prefilter,
        params_transigner_em,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

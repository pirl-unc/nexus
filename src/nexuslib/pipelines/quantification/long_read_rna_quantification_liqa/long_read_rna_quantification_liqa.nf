#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runLiqaQuantify } from '../../modules/liqa'

// Step 2. Input parameters
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.liqa_refgene_file = ''
// Optional arguments
params.params_liqa_quantify = '-max_distance 20 -f_weight 1'
params.delete_work_dir = false

if (params.params_liqa_quantify == true) {
    params_liqa_quantify = ''
} else {
    params_liqa_quantify = params.params_liqa_quantify
}

// Step 3. Print inputs and help
log.info """\
         ================================================
         Quantify RNA in long-read FASTQ files using LIQA
         ================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Quantify RNA in long-read FASTQ files using LIQA.

    usage: nexus run --nf-workflow long_read_rna_quantification_liqa.nf [required] [optional] [--help]

    required arguments:
        -c                                      :   Nextflow .config file.
        -w                                      :   Nextflow work directory path.
        --samples_tsv_file                      :   TSV file with the following columns:
                                                    'sample_id', 'bam_file', 'bam_bai_file'.
        --liqa_refgene_file                     :   LIQA refgene file.
        --output_dir                            :   Directory to which output files will be copied.

    optional arguments:
        --params_liqa_quantify                  :   LIQA quantify parameters (default: '"-max_distance 20 -f_weight 1"').
                                                    Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                       :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                        :   ${params.samples_tsv_file}
        liqa_refgene_file                       :   ${params.liqa_refgene_file}
        params_liqa_quantify                    :   ${params_liqa_quantify}
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
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow LONG_READ_RNA_QUANTIFICATION_LIQA {
    take:
        input_bam_files_ch
        liqa_refgene_file
        params_liqa_quantify
        output_dir
    main:
        runLiqaQuantify(
            input_bam_files_ch,
            liqa_refgene_file,
            params_liqa_quantify,
            output_dir
        )
    emit:
        runLiqaQuantify.out.f
}

workflow {
    LONG_READ_RNA_QUANTIFICATION_LIQA(
        input_bam_files_ch,
        params.liqa_refgene_file,
        params_liqa_quantify,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
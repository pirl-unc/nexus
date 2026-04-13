#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runLiqaIndex }        from '../../../tools/liqa'
include { runLiqaQuantify }     from '../../../tools/liqa'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir                   = ''
params.reference_genes_gtf_file     = ''

// Optional arguments
params.params_liqa_quantify         = '-max_distance 20 -f_weight 1'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow QUANTIFICATION_LIQA {
    take:
        input_bam_files_ch
        reference_genes_gtf_file
        params_liqa_quantify
        output_dir

    main:
        runLiqaIndex(reference_genes_gtf_file)

        runLiqaQuantify(
            input_bam_files_ch,
            runLiqaIndex.out.refgene_file,
            params_liqa_quantify,
            output_dir
        )

    emit:
        runLiqaQuantify.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ================================================
             Quantify RNA in long-read FASTQ files using LIQA
             ================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Quantify RNA in long-read FASTQ files using LIQA.

        usage: nexus run --nf-workflow quantification_liqa.nf [required] [optional] [--help]

        required arguments:
            -c                              :   Nextflow .config file.
            -w                              :   Nextflow work directory path.
            --samples_tsv_file              :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
            --reference_genes_gtf_file      :   Reference genes GTF file.
            --output_dir                    :   Directory to which output files will be copied.

        optional arguments:
            --params_liqa_quantify          :   LIQA quantify parameters (default: '"-max_distance 20 -f_weight 1"').
                                                Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_liqa_quantify = (params.params_liqa_quantify == true) ? '' : params.params_liqa_quantify

    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        reference_genes_gtf_file        :   ${params.reference_genes_gtf_file}
        params_liqa_quantify            :   ${params_liqa_quantify}
        output_dir                      :   ${params.output_dir}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    QUANTIFICATION_LIQA(
        input_bam_files_ch,
        params.reference_genes_gtf_file,
        params_liqa_quantify,
        params.output_dir
    )
}

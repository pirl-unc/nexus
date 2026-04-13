#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }            from '../../../tools/samtools'
include { runCircExplorer2Parse }       from '../../../tools/circexplorer2'
include { runCircExplorer2Annotate }    from '../../../tools/circexplorer2'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                                     = ''

// Required arguments
params.samples_tsv_file                         = ''
params.output_dir                               = ''
params.reference_genome_fasta_file              = ''
params.circexplorer2_gene_annotation_txt_file   = ''

// Optional arguments
params.params_circexplorer2_parse               = '-t STAR'
params.params_circexplorer2_annotate            = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_CIRCEXPLORER2 {
    take:
        input_files_ch             // channel: [val(sample_id), path(input_files_ch)]
        reference_genome_fasta_file
        circexplorer2_gene_annotation_txt_file
        params_circexplorer2_parse
        params_circexplorer2_annotate
        output_dir

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file      = runSamtoolsFaidx.out.fasta
        fasta_fai_file  = runSamtoolsFaidx.out.fai_file

        runCircExplorer2Parse(
            input_files_ch,
            params_circexplorer2_parse,
            output_dir
        )

        runCircExplorer2Annotate(
            runCircExplorer2Parse.out.f,
            fasta_file,
            fasta_fai_file,
            circexplorer2_gene_annotation_txt_file,
            params_circexplorer2_annotate,
            output_dir
        )

    emit:
        runCircExplorer2Annotate.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ==========================================
             Identify circular RNAs using CIRCexplorer2
             ==========================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run CIRCexplorer2 parse.
            2. Run CIRCexplorer2 annotate.

        usage: nexus run --nf-workflow variant_calling_circexplorer2.nf [required] [optional] [--help]

        required arguments:
            -c                                          :   Nextflow .config file.
            -w                                          :   Nextflow work directory path.
            --samples_tsv_file                          :   TSV file with the following columns:
                                                            'sample_id', 'input_file'.
            --output_dir                                :   Directory to which output files will be copied.
            --reference_genome_fasta_file               :   Reference genome FASTA file.
            --circexplorer2_gene_annotation_txt_file    :   CIRCexplorer2 gene annotation TXT file.

        optional arguments:
            --params_circexplorer2_parse                :   CIRCexplorer2 parse parameters (default: '"-t STAR"').
                                                            Note that the parameters need to be wrapped in quotes.
            --params_circexplorer2_annotate             :   CIRCexplorer2 annotate parameters (default: '""').
                                                            Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_circexplorer2_parse    = (params.params_circexplorer2_parse == true) ? '' : params.params_circexplorer2_parse
    def params_circexplorer2_annotate = (params.params_circexplorer2_annotate == true) ? '' : params.params_circexplorer2_annotate

    log.info"""\
        samples_tsv_file                            :   ${params.samples_tsv_file}
        output_dir                                  :   ${params.output_dir}
        reference_genome_fasta_file                 :   ${params.reference_genome_fasta_file}
        circexplorer2_gene_annotation_txt_file      :   ${params.circexplorer2_gene_annotation_txt_file}
        params_circexplorer2_parse                  :   ${params_circexplorer2_parse}
        params_circexplorer2_annotate               :   ${params.params_circexplorer2_annotate}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.input_file}") }
        .set { input_files_ch }

    VARIANT_CALLING_CIRCEXPLORER2(
        input_files_ch,
        params.reference_genome_fasta_file,
        params.circexplorer2_gene_annotation_txt_file,
        params_circexplorer2_parse,
        params_circexplorer2_annotate,
        params.output_dir
    )
}

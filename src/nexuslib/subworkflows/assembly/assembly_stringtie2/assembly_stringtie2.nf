#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runStringTie2 }                          from '../../../tools/stringtie2'
include { decompressFile as decompressGtf }        from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                     = ''

// Required arguments
params.samples_tsv_file         = ''
params.output_dir               = ''
params.reference_genes_gtf_file = ''

// Optional arguments
params.params_stringtie2    = '-L'

if (params.params_stringtie2 == true) {
    params_stringtie2 = ''
} else {
    params_stringtie2 = params.params_stringtie2
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         =================================================
         Assemble long-read RNA BAM files using StringTie2
         =================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Assemble transcripts using StringTie2.

    usage: nexus run --nf-workflow assembly_stringtie2.nf [required] [optional] [--help]

    required arguments:
        -c                          :   Nextflow .config file.
        -w                          :   Nextflow work directory path.
        --samples_tsv_file          :   TSV file with the following columns:
                                        'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                :   Directory to which output files will be copied.
        --reference_genes_gtf_file  :   Reference genes GTF file.

    optional arguments:
        --params_stringtie2         :   StringTie2 parameters (default: '"-L"').
                                        Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file            :   ${params.samples_tsv_file}
        output_dir                  :   ${params.output_dir}
        reference_genes_gtf_file    :   ${params.reference_genes_gtf_file}
        params_stringtie2           :   ${params_stringtie2}
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
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow ASSEMBLY_STRINGTIE2 {
    take:
        input_bam_files_ch            // channel: [val(sample_id), path(bam_file)]
        reference_genes_gtf_file
        params_stringtie2
        output_dir

    main:
        decompressGtf(reference_genes_gtf_file)

        runStringTie2(
            input_bam_files_ch,
            decompressGtf.out.f,
            params_stringtie2,
            output_dir
        )

    emit:
        runStringTie2.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ASSEMBLY_STRINGTIE2(
        input_bam_files_ch,
        params.reference_genes_gtf_file,
        params_stringtie2,
        params.output_dir
    )
}

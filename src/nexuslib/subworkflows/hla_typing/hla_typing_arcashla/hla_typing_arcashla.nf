#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runArcasHlaPairedEndMode }    from '../../../tools/arcashla'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow HLA_TYPING_ARCASHLA {
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

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ============================================================================
             Profile HLA alleles using paired-end RNA sequencing BAM files using arcasHLA
             ============================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Profile HLA alleles using paired-end RNA sequencing BAM files using arcasHLA.

        usage: nexus run --nf-workflow hla_typing_arcashla.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'bam_file', 'bam_bai_file'.
            --output_dir            :   Directory to which output files will be copied.
        """.stripIndent()
        exit 0
    }

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    HLA_TYPING_ARCASHLA(
        input_bam_files_ch,
        params.output_dir
    )
}

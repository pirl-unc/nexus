//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runRmatsBamMode }                    from '../../../tools/rmats'
include { decompressFile as decompressGtf }    from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir                   = ''
params.reference_genes_gtf_file     = ''

// Optional arguments
params.params_rmats                 = '-t paired --libType fr-unstranded --readLength 151'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ISOFORM_CHARACTERIZATION_RMATS {
    take:
        input_bam_files_ch
        reference_genes_gtf_file
        params_rmats
        output_dir

    main:
        decompressGtf(reference_genes_gtf_file)

        runRmatsBamMode(
            input_bam_files_ch,
            decompressGtf.out.f,
            params_rmats,
            output_dir
        )

    emit:
        runRmatsBamMode.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =================================
             Characterize isoforms using rMATS
             =================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run rMATS.

        usage: nexus run --nf-workflow isoform_characterization_rmats.nf [required] [optional] [--help]

        required arguments:
            -c                             :   Nextflow .config file.
            -w                             :   Nextflow work directory path.
            --samples_tsv_file             :   TSV file with the following columns:
                                               'sample_id', 'bam_file', 'bam_bai_file'.
            --output_dir                   :   Directory to which output files will be copied.
            --reference_genes_gtf_file     :   Reference genes GTF file.

        optional arguments:
            --params_rmats                 :   rMATS parameters (default: '"-t paired --libType fr-unstranded --readLength 151"').
                                               Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_rmats = (params.params_rmats == true) ? '' : params.params_rmats

    log.info"""\
        output_dir                     :   ${params.output_dir}
        reference_genes_gtf_file       :   ${params.reference_genes_gtf_file}
        params_rmats                   :   ${params_rmats}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    ISOFORM_CHARACTERIZATION_RMATS(
        input_bam_files_ch,
        params.reference_genes_gtf_file,
        params_rmats,
        params.output_dir
    )
}

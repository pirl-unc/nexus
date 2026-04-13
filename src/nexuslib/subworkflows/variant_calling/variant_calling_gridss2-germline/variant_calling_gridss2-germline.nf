#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runGridss2Germline }                     from '../../../tools/gridss2'
include { decompressFile as decompressFasta }      from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''

// Optional arguments
params.params_gridss2                   = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_GRIDSS2_GERMLINE {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        output_dir
        reference_genome_fasta_file
        params_gridss2

    main:
        decompressFasta(reference_genome_fasta_file)

        runGridss2Germline(
            input_bam_files_ch,
            decompressFasta.out.f,
            params_gridss2,
            output_dir
        )

    emit:
        runGridss2Germline.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ====================================================================================
             Identify germline variants in paired-end read DNA sequencing BAM files using GRIDSS2
             ====================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run GRIDSS2 (germline mode).

        usage: nexus run --nf-workflow variant_calling_gridss2-germline.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id',
                                                    'bam_file',
                                                    'bam_bai_file',
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.

        optional arguments:
            --params_gridss2                    :   GRIDSS2 gridss parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_gridss2 = (params.params_gridss2 == true) ? '' : params.params_gridss2

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_gridss2                      :   ${params_gridss2}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    VARIANT_CALLING_GRIDSS2_GERMLINE(
        input_bam_files_ch,
        params.output_dir,
        params.reference_genome_fasta_file,
        params_gridss2
    )
}

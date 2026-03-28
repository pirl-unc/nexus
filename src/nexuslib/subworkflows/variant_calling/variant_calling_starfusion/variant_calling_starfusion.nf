#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runStarFusion }       from '../../../tools/starfusion'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                     = ''

// Required arguments
params.samples_tsv_file         = ''
params.output_dir               = ''
params.genome_lib_dir           = ''

// Optional arguments
params.params_starfusion        = ''

if (params.params_starfusion == true) {
    params_starfusion = ''
} else {
    params_starfusion = params.params_starfusion
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         =======================================================================
         Identify fusion genes in STAR Chimeric junction files using STAR-Fusion
         =======================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run STAR-Fusion.

    usage: nexus run --nf-workflow variant_calling_starfusion.nf [required] [optional] [--help]

    required arguments:
        -c                      :   Nextflow .config file.
        -w                      :   Nextflow work directory path.
        --samples_tsv_file      :   TSV file with the following columns:
                                    'sample_id', 'fastq_file_1', 'fastq_file_2'.
        --output_dir            :   Directory to which output files will be copied.
        --genome_lib_dir        :   Genome lib directory.

    optional arguments:
        --params_starfusion     :   STAR-Fusion parameters (default: '""').
                                    Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        genome_lib_dir          :   ${params.genome_lib_dir}
        params_starfusion       :   ${params_starfusion}
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
        "${row.fastq_file_1}",
        "${row.fastq_file_2}") }
    .set { input_fastq_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_STAR_FUSION {
    take:
        input_fastq_files_ch              // channel: [val(sample_id), path(fastq_file_1), path(fastq_file_2)]
        genome_lib_dir
        params_starfusion
        output_dir

    main:
        runStarFusion(
            input_fastq_files_ch,
            genome_lib_dir,
            params_starfusion,
            output_dir
        )

    emit:
        runStarFusion.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_STAR_FUSION(
        input_fastq_files_ch,
        params.genome_lib_dir,
        params_starfusion,
        params.output_dir
    )
}

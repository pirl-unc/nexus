#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runNexusFilterRNABloom2Transcripts }  from '../../../tools/nexus'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_filter        = '--min-mapping-quality 30 --min-read-support 3 --min-fraction-match 0.5'

if (params.params_filter == true) {
    params_filter = ''
} else {
    params_filter = params.params_filter
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ========================================
         Filter RNABloom2 transcripts using Nexus
         ========================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run filter_rnabloom2_transcripts (Nexus).

    usage: nexus run --nf-workflow utilities_filter-rnabloom2-transcripts.nf [required] [optional] [--help]

    required arguments:
        -c                      :   Nextflow .config file.
        -w                      :   Nextflow work directory path.
        --samples_tsv_file      :   TSV file with the following columns: 'sample_id', 'rnabloom2_output_dir''.
        --output_dir            :   Directory to which output files will be copied.

    optional arguments:
        --params_filter         :   filter_rnabloom2_transcripts parameters (default: '"--min-mapping-quality 30 --min-read-support 3 --min-fraction-match 0.5"').
                                    Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_filter           :   ${params_filter}
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
        "${row.rnabloom2_output_dir}") }
    .set { input_dir_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow UTILITIES_FILTER_RNABLOOM2_TRANSCRIPTS {
    take:
        input_dir_ch             // channel: [val(sample_id), path(rnabloom2_output_dir)]
        params_filter
        output_dir

    main:
        runNexusFilterRNABloom2Transcripts(
            input_dir_ch,
            params_filter,
            output_dir
        )

    emit:
        runNexusFilterRNABloom2Transcripts.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    UTILITIES_FILTER_RNABLOOM2_TRANSCRIPTS(
        input_dir_ch,
        params_filter,
        params.output_dir
    )
}

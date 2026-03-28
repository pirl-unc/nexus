#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runHifiasm }      from '../../../tools/hifiasm'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_hifiasm       = ''

if (params.params_hifiasm == true) {
    params_hifiasm = ''
} else {
    params_hifiasm = params.params_hifiasm
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ==============================================
         Assemble long-read DNA BAM files using hifiasm
         ==============================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Perform de novo assembly using hifiasm.

    usage: nexus run --nf-workflow assembly_hifiasm.nf [required] [optional] [--help]

    required arguments:
        -c                      :   Nextflow .config file.
        -w                      :   Nextflow work directory path.
        --samples_tsv_file      :   TSV file with the following columns:
                                    'sample_id', 'fastq_file_1', 'fastq_file_2', ..., 'fastq_file_*'.
        --output_dir            :   Directory to which output files will be copied.

    optional arguments:
        --params_hifiasm        :   hifiasm parameters (default: '""').
                                    Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_hifiasm          :   ${params_hifiasm}
    """.stripIndent()
}

// ------------------------------------------------------------
// Step 4. Set channels
// ------------------------------------------------------------
Channel
    .fromPath(params.samples_tsv_file)
    .splitCsv(header: true, sep: '\t')
    .map { row ->
        def sample_id = row.sample_id

        // Find all columns whose names start with "fastq_file"
        def fastq_cols = row.keySet().findAll { it.startsWith('fastq_file') }

        // Collect non-empty paths as files
        def fastqs = fastq_cols
            .collect { row[it] }
            .findAll { it }
            .collect { file(it) }

        tuple(sample_id, fastqs)
    }
    .set { input_fastq_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow ASSEMBLY_HIFIASM {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), [path(fastq_file_1), path(fastq_file_2), ...,]]
        params_hifiasm
        output_dir

    main:
        runHifiasm(
            input_fastq_files_ch,
            params_hifiasm,
            output_dir
        )

    emit:
        runHifiasm.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ASSEMBLY_HIFIASM(
        input_fastq_files_ch,
        params_hifiasm,
        params.output_dir
    )
}

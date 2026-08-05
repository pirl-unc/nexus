#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSpecImmune }    from '../../../tools/specimmune'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_specimmune    = '-i HLA -y pacbio-hifi --seq_tech rna --RNA_type traditional'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow HLA_TYPING_SPECIMMUNE {
    take:
        input_fastq_files_ch
        params_specimmune
        output_dir

    main:
        runSpecImmune(
            input_fastq_files_ch,
            params_specimmune,
            output_dir
        )
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ==============================================================================================
             Profile HLA / KIR / CYP / IG-TR alleles from long-read DNA or RNA FASTQ files using SpecImmune
             ==============================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Profile HLA / KIR / CYP / IG_TR / extend alleles from long-read
               FASTQ files using SpecImmune (scripts/main.py).

        usage: nexus run --nf-workflow hla_typing_specimmune.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_specimmune     :   SpecImmune (scripts/main.py) extra CLI parameters
                                        (default: '"-i HLA -y pacbio"').
                                        Use -i to select the typing target
                                        (HLA | KIR | CYP | IG_TR | extend).
                                        Use -y to select the read type
                                        (nanopore | pacbio | pacbio-hifi).
                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_specimmune = (params.params_specimmune == true) ? '' : params.params_specimmune

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_specimmune       :   ${params_specimmune}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_fastq_files_ch }

    HLA_TYPING_SPECIMMUNE(
        input_fastq_files_ch,
        params_specimmune,
        params.output_dir
    )
}

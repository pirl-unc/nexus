#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSalmonIndex }          from '../../../tools/salmon'
include { runSalmonFastqMode }      from '../../../tools/salmon'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                                 = ''

// Required arguments
params.samples_tsv_file                     = ''
params.output_dir                           = ''
params.reference_transcripts_fasta_file     = ''

// Optional arguments
params.params_salmon_index                  = '--gencode'
params.params_salmon_quant                  = '--libType IU --seqBias --gcBias --posBias'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow QUANTIFICATION_SALMON_FASTQ {
    take:
        input_fastq_files_ch
        reference_transcripts_fasta_file
        params_salmon_index
        params_salmon_quant
        output_dir

    main:
        runSalmonIndex(
            reference_transcripts_fasta_file,
            params_salmon_index
        )

        runSalmonFastqMode(
            input_fastq_files_ch,
            reference_transcripts_fasta_file,
            runSalmonIndex.out.f,
            params_salmon_quant,
            output_dir
        )

    emit:
        runSalmonFastqMode.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===================================================
             Quantify RNA in paired-end FASTQ files using Salmon
             ===================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Quantify RNA in paired-end FASTQ files using Salmon.

        usage: nexus run --nf-workflow quantification_salmon-fastq.nf [required] [optional] [--help]

        required arguments:
            -c                                      :   Nextflow .config file.
            -w                                      :   Nextflow work directory path.
            --samples_tsv_file                      :   TSV file with the following columns:
                                                        'sample_id', 'fastq_file_1', 'fastq_file_2'.
            --output_dir                            :   Directory to which output files will be copied.
            --reference_transcripts_fasta_file      :   Reference transcripts FASTA file.

        optional arguments:
            --params_salmon_index                   :   Salmon index parameters (default: '"--gencode "').
                                                        Note that the parameters need to be wrapped in quotes.
            --params_salmon_quant                   :   Salmon parameters (default: '"--libType IU --seqBias --gcBias --posBias"').
                                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_salmon_index = (params.params_salmon_index == true) ? '' : params.params_salmon_index
    def params_salmon_quant = (params.params_salmon_quant == true) ? '' : params.params_salmon_quant

    log.info"""\
        samples_tsv_file                        :   ${params.samples_tsv_file}
        output_dir                              :   ${params.output_dir}
        reference_transcripts_fasta_file        :   ${params.reference_transcripts_fasta_file}
        params_salmon_index                     :   ${params_salmon_index}
        params_salmon_quant                     :   ${params_salmon_quant}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file_1}",
            "${row.fastq_file_2}") }
        .set { input_fastq_files_ch }

    QUANTIFICATION_SALMON_FASTQ(
        input_fastq_files_ch,
        params.reference_transcripts_fasta_file,
        params_salmon_index,
        params_salmon_quant,
        params.output_dir
    )
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runKallistoIndex }                    from '../../../tools/kallisto'
include { runKallistoQuantPairedEndReads }      from '../../../tools/kallisto'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                                 = ''

// Required arguments
params.samples_tsv_file                     = ''
params.output_dir                           = ''
params.reference_transcripts_fasta_file     = ''

// Optional arguments
params.params_kallisto_index                = '-k 31'
params.params_kallisto_quant = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow QUANTIFICATION_KALLISTO_PE {
    take:
        input_fastq_files_ch
        reference_transcripts_fasta_file
        params_kallisto_quant
        output_dir

    main:
        runKallistoIndex(
            reference_transcripts_fasta_file,
            params_kallisto_index
        )

        runKallistoQuantPairedEndReads(
            input_fastq_files_ch,
            runKallistoIndex.out.index_file,
            params_kallisto_quant,
            output_dir
        )

    emit:
        runKallistoQuantPairedEndReads.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===========================================================
             Quantify RNA in pairend-end read FASTQ files using Kallisto
             ===========================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Quantify RNA in paired-end read FASTQ files using Kallisto.

        usage: nexus run --nf-workflow quantification_kallisto-pe.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'fastq_file_1', 'fastq_file_2'.
            --reference_transcripts_fasta_file  :   Reference transcripts FASTA file.
            --output_dir                        :   Directory to which output files will be copied.

        optional arguments:
            --params_kallisto_quant             :   Kallisto quant parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_kallisto_index = (params.params_kallisto_index == true) ? '' : params.params_kallisto_index
    def params_kallisto_quant = (params.params_kallisto_quant == true) ? '' : params.params_kallisto_quant

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        reference_transcripts_fasta_file    :   ${params.reference_transcripts_fasta_file}
        params_kallisto_quant               :   ${params_kallisto_quant}
        output_dir                          :   ${params.output_dir}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file_1}",
            "${row.fastq_file_2}") }
        .set { input_fastq_files_ch }

    QUANTIFICATION_KALLISTO_PE(
        input_fastq_files_ch,
        params.reference_transcripts_fasta_file,
        params_kallisto_quant,
        params.output_dir
    )
}

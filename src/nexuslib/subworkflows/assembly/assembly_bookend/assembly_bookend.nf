#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }               from '../../../tools/samtools'
include { decompressFile as decompressFasta }   from '../../../tools/utils'
include { runBookendAssemble }                  from '../../../tools/bookend'
include { runBookendFasta }                     from '../../../tools/bookend'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                     = ''

// Required arguments
params.samples_tsv_file         = ''
params.output_dir               = ''

// Optional arguments
params.params_bookend_assemble  = ''
params.params_bookend_fasta     = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ASSEMBLY_BOOKEND {
    take:
        input_bam_files_ch            // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        params_bookend_assemble
        params_bookend_fasta
        output_dir

    main:
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidxFasta(decompressFasta.out.f)
        uncompressed_fasta     = runSamtoolsFaidxFasta.out.fasta
        uncompressed_fasta_fai = runSamtoolsFaidxFasta.out.fai_file

        runBookendAssemble(
            input_bam_files_ch,
            params_bookend_assemble,
            output_dir
        )

        runBookendFasta(
            runBookendAssemble.out.f,
            uncompressed_fasta,
            uncompressed_fasta_fai,
            params_bookend_fasta,
            output_dir
        )

    emit:
        runBookendFasta.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ================================================
             Assemble long-read RNA FASTQ files using Bookend
             ================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Assemble transcripts using Bookend.

        usage: nexus run --nf-workflow assembly_bookend.nf [required] [optional] [--help]

        required arguments:
            -c                              :   Nextflow .config file.
            -w                              :   Nextflow work directory path.
            --samples_tsv_file              :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
            --output_dir                    :   Directory to which output files will be copied.
            --reference_genome_fasta_file   :   Reference genome FASTA file.

        optional arguments:
            --params_bookend_assemble       :   bookend assemble parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
            --params_bookend_fasta          :   bookend fasta parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_bookend_assemble    = (params.params_bookend_assemble == true) ? '' : params.params_bookend_assemble
    def params_bookend_fasta       = (params.params_bookend_fasta == true) ? '' : params.params_bookend_fasta

    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        reference_genome_fasta_file     :   ${params.reference_genome_fasta_file}
        params_bookend_assemble         :   ${params_bookend_assemble}
        params_bookend_fasta            :   ${params_bookend_fasta}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    ASSEMBLY_BOOKEND(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params_bookend_assemble,
        params_bookend_fasta,
        params.output_dir
    )
}

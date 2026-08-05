#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { decompressFile as decompressFasta } from '../../../tools/utils'
include { runSamtoolsFaidxFasta }             from '../../../tools/samtools'
include { runGatk4CreateSequenceDictionary }  from '../../../tools/gatk4'
include { runIsolaser }                       from '../../../tools/isolaser'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''

// Optional arguments
params.params_isolaser                  = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_ISOLASER {
    take:
        input_files_ch              // channel: [val(sample_id), path(bam_file), path(bam_bai_file), path(fastq_file), path(gtf_file)]
        reference_genome_fasta_file
        params_isolaser
        output_dir

    main:
        // Decompress + index the reference FASTA. isoLASER (pyfaidx) needs an
        // UNCOMPRESSED .fa, so do not bgzip it.
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidxFasta(decompressFasta.out.f)
        fasta_file      = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file  = runSamtoolsFaidxFasta.out.fai_file

        runGatk4CreateSequenceDictionary(fasta_file)
        fasta_dict_file = runGatk4CreateSequenceDictionary.out.f

        runIsolaser(
            input_files_ch,
            fasta_file,
            fasta_fai_file,
            fasta_dict_file,
            params_isolaser,
            output_dir
        )

    emit:
        runIsolaser.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ==========================================================================
             Identify RNA variants in long-read RNA sequencing BAM files using Isolaser
             ==========================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run Isolaser.

        usage: nexus run --nf-workflow variant_calling_isolaser.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'bam_file', 'bam_bai_file', 'fastq_file', 'gtf_file' (from Talon, Bambu or ESPRESSO).
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.

        optional arguments:
            --params_isolaser                   :   Isolaser parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_isolaser = (params.params_isolaser == true) ? '' : params.params_isolaser

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_isolaser                     :   ${params_isolaser}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}",
            "${row.fastq_file}",
            "${row.gtf_file}") }
        .set { input_files_ch }

    VARIANT_CALLING_ISOLASER(
        input_files_ch,
        params.reference_genome_fasta_file,
        params_isolaser,
        params.output_dir
    )
}


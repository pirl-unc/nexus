#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { decompressFile as decompressFasta }       from '../../../tools/utils'
include { runSamtoolsFaidxFasta }                   from '../../../tools/samtools'
include { runFlairAlign }                           from '../../../tools/flair'
include { runLongshot }                             from '../../../tools/longshot'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''

// Optional arguments
params.params_flair_align               = ''
params.params_longshot                  = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow HAPLOTAGGING_FLAIR_LONGSHOT {
    take:
        input_fastq_files_ch             // channel: [val(sample_id), path(fastq_file)]
        reference_genome_fasta_file
        params_flair_align
        params_longshot
        output_dir

    main:
        // Step 1. Decompress and index reference genome FASTA file
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidxFasta(decompressFasta.out.f)
        fasta_file          = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file      = runSamtoolsFaidxFasta.out.fai_file

        // Step 2. Run Flair
        runFlairAlign(
            input_fastq_files_ch,
            fasta_file,
            params_flair_align,
            output_dir
        )

        // Step 3. Run LongShot
        // runFlairAlign emits (sample_id, bam, bam.bai, bed); drop the bed for longshot input
        longshot_input_ch = runFlairAlign.out.f
            .map { sample_id, bam_file, bam_bai_file, bed_file -> tuple(sample_id, bam_file, bam_bai_file) }

        runLongshot(
            longshot_input_ch,
            fasta_file,
            fasta_fai_file,
            params_longshot,
            output_dir
        )

    emit:
        // channel: [val(sample_id), path(vcf_file), path(bam_file), path(bam_bai_file)]
        runLongshot.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =================================================================
             Haplotype long-read DNA sequencing BAM files using Flair+Longshot
             =================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run flair align.
            2. Run longshot.

        usage: nexus run --nf-workflow haplotagging_flair-longshot.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id',
                                                    'fastq_file'
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.

        optional arguments:
            --params_flair_align                :   flair align parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
            --params_longshot                   :   longshot parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_flair_align = (params.params_flair_align == true) ? '' : params.params_flair_align
    def params_longshot = (params.params_longshot == true) ? '' : params.params_longshot

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        preset                              :   ${params.preset}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_flair_align                  :   ${params_flair_align}
        params_longshot                     :   ${params_longshot}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_fastq_files_ch }

    HAPLOTAGGING_FLAIR_LONGSHOT(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params_flair_align,
        params_longshot,
        params.output_dir
    )
}

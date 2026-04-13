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
include { runLongcallR }                            from '../../../tools/longcallr'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.reference_genes_gtf_file         = ''
params.preset                           = ''

// Optional arguments
params.params_longcallr                 = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow HAPLOTYPING_LONGCALLR {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        preset
        reference_genome_fasta_file
        reference_genes_gtf_file
        params_longcallr
        output_dir

    main:
        // Step 1. Decompress and index reference genome FASTA file
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidxFasta(decompressFasta.out.f)
        fasta_file          = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file      = runSamtoolsFaidxFasta.out.fai_file

        // Step 2. Run LongcallR
        runLongcallR(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            reference_genes_gtf_file,
            preset,
            params_longcallr,
            output_dir
        )

    emit:
        runLongcallR.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ============================================================
             Haplotype long-read DNA sequencing BAM files using LongcallR
             ============================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run 'longcallr' command.

        usage: nexus run --nf-workflow haplotyping_longcallr.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id',
                                                    'bam_file',
                                                    'bam_bai_file'
            --output_dir                        :   Directory to which output files will be copied.
            --preset                            :   LongcallR preset (choices: hifi-isoseq, hifi-masseq, ont-cdna, ont-drna).
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --reference_genes_gtf_file          :   Reference genes GTF file.

        optional arguments:
            --params_longcallr                  :   longcallr parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_longcallr = (params.params_longcallr == true) ? '' : params.params_longcallr

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        preset                              :   ${params.preset}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_longcallr                    :   ${params_longcallr}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    HAPLOTYPING_LONGCALLR(
        input_bam_files_ch,
        params.preset,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params_longcallr,
        params.output_dir
    )
}

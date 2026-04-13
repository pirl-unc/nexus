#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }            from '../../../tools/samtools'
include { runBambuDiscoverQuantify }    from '../../../tools/bambu'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.reference_genes_gtf_file         = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow QUANTIFICATION_BAMBU {
    take:
        input_bam_files_ch
        reference_genome_fasta_file
        reference_genes_gtf_file
        output_dir

    main:
        // Step 1. Index reference genome FASTA file
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file          = runSamtoolsFaidx.out.fasta
        fasta_fai_file      = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file      = runSamtoolsFaidx.out.gzi_file

        // Step 2. Run Bambu
        runBambuDiscoverQuantify(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            reference_genes_gtf_file,
            output_dir
        )

    emit:
        runBambuDiscoverQuantify.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =================================================
             Quantify RNA in long-read FASTQ files using Bambu
             =================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Quantify RNA in long-read FASTQ files using Bambu.

        usage: nexus run --nf-workflow quantification_bambu.nf [required] [optional] [--help]

        required arguments:
            -c                              :   Nextflow .config file.
            -w                              :   Nextflow work directory path.
            --samples_tsv_file              :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
            --output_dir                    :   Directory to which output files will be copied.
            --reference_genome_fasta_file   :   Reference genome FASTA file.
            --reference_genes_gtf_file      :   Reference genes GTF file.
        """.stripIndent()
        exit 0
    }

    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        reference_genome_fasta_file     :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file        :   ${params.reference_genes_gtf_file}
        output_dir                      :   ${params.output_dir}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    QUANTIFICATION_BAMBU(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params.output_dir
    )
}

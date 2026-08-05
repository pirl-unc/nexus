//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runFlair3Transcriptome }             from '../../../tools/flair3'
include { decompressFile as decompressFasta }  from '../../../tools/utils'
include { decompressFile as decompressGtf }    from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir                   = ''
params.reference_genome_fasta_file  = ''
params.reference_genes_gtf_file     = ''

// Optional arguments
params.params_flair_transcriptome   = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ISOFORM_CHARACTERIZATION_FLAIR3 {
    take:
        input_bam_files_ch            // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        reference_genes_gtf_file
        params_flair_transcriptome
        output_dir

    main:
        // Step 1. Decompress reference files if needed
        decompressFasta(reference_genome_fasta_file)
        decompressGtf(reference_genes_gtf_file)

        // Step 2. Run Flair align
        runFlair3Transcriptome(
            input_bam_files_ch,
            decompressFasta.out.f,
            decompressGtf.out.f,
            params_flair_transcriptome,
            output_dir
        )

    emit:
        runFlair3Transcriptome.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ==================================
             Characterize isoforms using Flair3
             ==================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run flair transcriptome.

        usage: nexus run --nf-workflow isoform_characterization_flair3.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'bam_file', 'bam_bai_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --reference_genes_gtf_file          :   Reference genes GTF file.

        optional arguments:
            --params_flair_transcriptome        :   Flair3 'transcriptome' parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_flair_transcriptome    = (params.params_flair_transcriptome == true) ? '' : params.params_flair_transcriptome

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_flair_transcriptome          :   ${params_flair_transcriptome}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    ISOFORM_CHARACTERIZATION_FLAIR3(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params_flair_transcriptome,
        params.output_dir
    )
}

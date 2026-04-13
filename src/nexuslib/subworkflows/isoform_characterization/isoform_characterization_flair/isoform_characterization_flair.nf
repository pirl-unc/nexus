//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runFlairAlign }                      from '../../../tools/flair'
include { runFlairCorrect }                    from '../../../tools/flair'
include { runFlairCollapse }                   from '../../../tools/flair'
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
params.params_flair_align           = ''
params.params_flair_correct         = ''
params.params_flair_collapse        = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ISOFORM_CHARACTERIZATION_FLAIR {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        reference_genome_fasta_file
        reference_genes_gtf_file
        params_flair_align
        params_flair_correct
        params_flair_collapse
        output_dir

    main:
        // Step 1. Decompress reference files if needed
        decompressFasta(reference_genome_fasta_file)
        decompressGtf(reference_genes_gtf_file)

        // Step 2. Run Flair align
        runFlairAlign(
            input_fastq_files_ch,
            decompressFasta.out.f,
            params_flair_align,
            output_dir
        )

        // Step 3. Run Flair correct
        runFlairCorrect(
            runFlairAlign.out.f,
            decompressFasta.out.f,
            decompressGtf.out.f,
            params_flair_correct,
            output_dir
        )

        // Step 4. Run Flair collapse
        runFlairCorrect.out.f.set{ run_flair_correct_output_ch }
        run_flair_correct_output_ch
           .join(input_fastq_files_ch)
           .set{ run_flair_collapse_input_ch }
        runFlairCollapse(
            run_flair_collapse_input_ch,
            decompressFasta.out.f,
            decompressGtf.out.f,
            params_flair_collapse,
            output_dir
        )

    emit:
        runFlairCollapse.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =================================
             Characterize isoforms using Flair
             =================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run flair align.
            2. Run flair correct.
            3. Run flair collapse.

        usage: nexus run --nf-workflow isoform_characterization_flair.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'fastq_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --reference_genes_gtf_file          :   Reference genes GTF file.

        optional arguments:
            --params_flair_align                :   Flair 'align' parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
            --params_flair_correct              :   Flair 'correct' parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
            --params_flair_collapse             :   Flair 'collapse' parameters (default: '" "').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_flair_align    = (params.params_flair_align == true) ? '' : params.params_flair_align
    def params_flair_correct  = (params.params_flair_correct == true) ? '' : params.params_flair_correct
    def params_flair_collapse = (params.params_flair_collapse == true) ? '' : params.params_flair_collapse

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_flair_align                  :   ${params_flair_align}
        params_flair_correct                :   ${params_flair_correct}
        params_flair_collapse               :   ${params_flair_collapse}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_fastq_files_ch }

    ISOFORM_CHARACTERIZATION_FLAIR(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params_flair_align,
        params_flair_correct,
        params_flair_collapse,
        params.output_dir
    )
}

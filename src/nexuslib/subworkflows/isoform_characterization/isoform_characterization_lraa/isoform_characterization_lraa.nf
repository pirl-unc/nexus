//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }    from '../../../tools/samtools'
include { runLRAA }             from '../../../tools/lraa'

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
params.params_lraa                  = '--HiFi'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ISOFORM_CHARACTERIZATION_LRAA {
    take:
        input_bam_files_ch
        reference_genome_fasta_file
        reference_genes_gtf_file
        params_lraa
        output_dir

    main:
        // Step 1. Index reference genome FASTA file
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file          = runSamtoolsFaidx.out.fasta
        fasta_fai_file      = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file      = runSamtoolsFaidx.out.gzi_file

        // Step 2. Run LRAA
        runLRAA(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            reference_genes_gtf_file,
            params_lraa,
            output_dir
        )

    emit:
        runLRAA.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ================================
             Characterize isoforms using LRAA
             ================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run LRAA.

        usage: nexus run --nf-workflow isoform_characterization_lraa.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'bam_file', 'bam_bai_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --reference_genes_gtf_file          :   Reference genes GTF file.

        optional arguments:
            --params_lraa                       :   LRAA parameters (default: '"--HiFi"').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_lraa = (params.params_lraa == true) ? '' : params.params_lraa

    log.info"""\
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_lraa                         :   ${params_lraa}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    ISOFORM_CHARACTERIZATION_LRAA(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params_lraa,
        params.output_dir
    )
}

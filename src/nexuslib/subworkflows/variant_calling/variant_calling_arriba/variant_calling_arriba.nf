#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }    from '../../../tools/samtools'
include { runArriba }           from '../../../tools/arriba'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.reference_genes_gtf_file         = ''
params.protein_domains_gff3_file        = ''

// Optional arguments
params.params_arriba                    = '-S 3 -f blacklist -i chr*'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_ARRIBA {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        reference_genes_gtf_file
        protein_domains_gff3_file
        params_arriba
        output_dir

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file          = runSamtoolsFaidx.out.fasta
        fasta_fai_file      = runSamtoolsFaidx.out.fai_file

        runArriba(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            reference_genes_gtf_file,
            protein_domains_gff3_file,
            params_arriba,
            output_dir
        )

    emit:
        runArriba.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ==========================================================================
             Identify fusion genes in paired-read RNA sequencing BAM files using Arriba
             ==========================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Align paired-end reads to a reference genome index using STAR.
            2. Run Arriba.

        usage: nexus run --nf-workflow variant_calling_arriba.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'bam_file', 'bam_bai_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --reference_genes_gtf_file          :   Reference genes GTF file.
            --protein_domains_gff3_file         :   Protein domains GFF3 file.

        optional arguments:
            --params_arriba                     :   Arriba parameters (default: '"-S 3 -f blacklist -i chr*"').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_arriba = (params.params_arriba == true) ? '' : params.params_arriba

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        protein_domains_gff3_file           :   ${params.protein_domains_gff3_file}
        params_arriba                       :   ${params_arriba}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    VARIANT_CALLING_ARRIBA(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params.protein_domains_gff3_file,
        params_arriba,
        params.output_dir
    )
}

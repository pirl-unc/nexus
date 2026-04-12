#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }    from '../../../tools/samtools'
include { runPbfusion }         from '../../../tools/pbfusion'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.reference_genes_gtf_file         = ''

// Optional arguments
params.params_pbfusion_discover         = '--min-coverage 3 --min-mean-mapq 20 --gtf-transcript-allow-lncRNA'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_PBFUSION {
    take:
        input_bam_files_ch              // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        reference_genes_gtf_file
        params_pbfusion_discover
        output_dir

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file                  = runSamtoolsFaidx.out.fasta
        fasta_fai_file              = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file              = runSamtoolsFaidx.out.gzi_file

        runPbfusion(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            reference_genes_gtf_file,
            params_pbfusion_discover,
            output_dir
        )

    emit:
        runPbfusion.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ==========================================================================
             Identify RNA variants in long-read RNA sequencing BAM files using Pbfusion
             ==========================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run Pbfusion.

        usage: nexus run --nf-workflow variant_calling_pbfusion.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id',
                                                    'bam_file',
                                                    'bam_bai_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --reference_genes_gtf_file          :   Reference genes GTF file.

        optional arguments:
            --params_pbfusion_discover          :   Pbfusion discover parameters (default: '"--min-coverage 3 --min-mean-mapq 20 --gtf-transcript-allow-lncRNA"').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_pbfusion_discover = (params.params_pbfusion_discover == true) ? '' : params.params_pbfusion_discover

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_pbfusion_discover            :   ${params_pbfusion_discover}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    VARIANT_CALLING_PBFUSION(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params_pbfusion_discover,
        params.output_dir
    )
}

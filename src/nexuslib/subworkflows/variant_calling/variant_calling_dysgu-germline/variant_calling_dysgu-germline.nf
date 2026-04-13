#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }            from '../../../tools/samtools'
include { runDysguGermlineMode }        from '../../../tools/dysgu'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''

// Optional arguments
params.params_dysgu_run                 = '--mode pacbio-revio --min-support 3 --min-size 30 --mq 20'
params.params_dysgu_filter              = '--support-fraction 0.05 --min-mapq 20 --pass-prob 0.2'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_DYSGU_GERMLINE {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        params_dysgu_run
        params_dysgu_filter
        output_dir

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file      = runSamtoolsFaidx.out.fasta
        fasta_fai_file  = runSamtoolsFaidx.out.fai_file

        runDysguGermlineMode(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            params_dysgu_run,
            params_dysgu_filter,
            output_dir
        )

    emit:
        runDysguGermlineMode.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ==========================================================================================================
             Identify germline structural variants in long-read or paired-end read DNA sequencing BAM files using Dysgu
             ==========================================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run Dysgu.

        usage: nexus run --nf-workflow variant_calling_dysgu-germline.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'bam_file', 'bam_bai_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.

        optional arguments:
            --params_dysgu_run                  :   Dysgu run parameters (default: '"--mode pacbio-revio --min-support 3 --min-size 30 --mq 20"').
                                                    Note that the parameters need to be wrapped in quotes.
            --params_dysgu_filter               :   Dysgu filter parameters (default: '"--support-fraction 0.05 --min-mapq 20 --pass-prob 0.2"').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_dysgu_run = (params.params_dysgu_run == true) ? '' : params.params_dysgu_run
    def params_dysgu_filter = (params.params_dysgu_filter == true) ? '' : params.params_dysgu_filter

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_dysgu_run                    :   ${params_dysgu_run}
        params_dysgu_filter                 :   ${params_dysgu_filter}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    VARIANT_CALLING_DYSGU_GERMLINE(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params_dysgu_run,
        params_dysgu_filter,
        params.output_dir
    )
}

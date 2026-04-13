#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }                from '../../../tools/samtools'
include { runDelly2ShortReadGermline }      from '../../../tools/delly2'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir                   = ''
params.reference_genome_fasta_file  = ''
params.exclude_tsv_file             = ''

// Optional arguments
params.params_delly2call            = '--map-qual 20'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_DELLY2_SHORTREAD_GERMLINE {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        output_dir
        reference_genome_fasta_file
        exclude_tsv_file
        params_delly2call

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file      = runSamtoolsFaidx.out.fasta
        fasta_fai_file  = runSamtoolsFaidx.out.fai_file

        runDelly2ShortReadGermline(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            exclude_tsv_file,
            params_delly2call,
            output_dir
        )

    emit:
        runDelly2ShortReadGermline.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===================================================================================
             Identify germline variants in paired-end read DNA sequencing BAM files using Delly2
             ===================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run Delly2 (germline mode).

        usage: nexus run --nf-workflow variant_calling_delly2-sr-germline.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id',
                                                    'bam_file',
                                                    'bam_bai_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --exclude_tsv_file                  :   Delly2 'call' --exclude TSV file.

        optional arguments:
            --params_delly2call                 :   Delly2 'call' parameters (default: '"--map-qual 20"').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_delly2call = (params.params_delly2call == true) ? '' : params.params_delly2call

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        exclude_tsv_file                    :   ${params.exclude_tsv_file}
        params_delly2call                   :   ${params_delly2call}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    VARIANT_CALLING_DELLY2_SHORTREAD_GERMLINE(
        input_bam_files_ch,
        params.output_dir,
        params.reference_genome_fasta_file,
        params.exclude_tsv_file,
        params_delly2call
    )
}

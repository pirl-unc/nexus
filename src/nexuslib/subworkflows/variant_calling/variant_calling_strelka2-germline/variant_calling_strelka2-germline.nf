#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }            from '../../../tools/samtools'
include { runStrelka2GermlineMode }     from '../../../tools/strelka2'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''

// Optional arguments
params.params_strelka2                  = ''

if (params.params_strelka2 == true) {
    params_strelka2 = ''
} else {
    params_strelka2 = params.params_strelka2
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         =====================================================================================
         Identify germline variants in paired-end read DNA sequencing BAM files using Strelka2
         =====================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1.  Run Strelka2.

    usage: nexus run --nf-workflow variant_calling_strelka2-germline.nf [required] [optional] [--help]

    required arguments:
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be symlinked.
        --reference_genome_fasta_file       :   Reference genome FASTA file.

    optional arguments:
        --params_strelka2                   :   Strelka2 parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_strelka2                     :   ${params_strelka2}
    """.stripIndent()
}

// ------------------------------------------------------------
// Step 4. Set channels
// ------------------------------------------------------------
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_STRELKA2_GERMLINE {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        params_strelka2
        output_dir

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file                  = runSamtoolsFaidx.out.fasta
        fasta_fai_file              = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file              = runSamtoolsFaidx.out.gzi_file

        runStrelka2GermlineMode(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            fasta_gzi_file,
            params_strelka2,
            output_dir
        )

    emit:
        runStrelka2GermlineMode.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_STRELKA2_GERMLINE(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params_strelka2,
        params.output_dir
    )
}

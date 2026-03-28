#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }    from '../../../tools/samtools'
include { runClair3 }                from '../../../tools/clair3'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''

// Optional arguments
params.params_clair3                    = '--model_path=/opt/models/hifi_sequel2/ --platform=hifi --min_coverage=3'

if (params.params_clair3 == true) {
    params_clair3 = ''
} else {
    params_clair3 = params.params_clair3
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ========================================================================================================================================
         Identify DNA variants in long or paired-end read DNA sequencing BAM files using Clair3 (Zheng et al., Nature Computational Science 2022)
         ========================================================================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run Clair3.

    usage: nexus run --nf-workflow variant_calling_clair3.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.

    optional arguments:
        --params_clair3                     :   Clair3 parameters (default: '"--model_path=/opt/models/hifi_sequel2/ --platform=hifi --min_coverage=3"').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_clair3                       :   ${params_clair3}
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
workflow VARIANT_CALLING_CLAIR3 {
    take:
        input_bam_files_ch              // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        params_clair3
        output_dir

    main:
        runSamtoolsFaidxFasta(reference_genome_fasta_file)
        fasta_file      = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file  = runSamtoolsFaidxFasta.out.fai_file

        runClair3(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            params_clair3,
            output_dir
        )

    emit:
        runClair3.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_CLAIR3(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params_clair3,
        params.output_dir
    )
}


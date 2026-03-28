#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }       from '../../../tools/samtools'
include { runHimutCall }                from '../../../tools/himut'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.region_list_file                 = ''

// Optional arguments
params.params_himut_call                = '--min_qv 20 --min_mapq 20 --min_bq 20'

if (params.params_himut_call == true) {
    params_himut_call = ''
} else {
    params_himut_call = params.params_himut_call
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ==================================================================================================
         Identify DNA variants in long-read DNA sequencing BAM files using Himut (Lee et al., bioRxiv 2025)
         ==================================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run himut call.

    usage: nexus run --nf-workflow variant_calling_himut.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.
        --region_list_file                  :   Region list file (chromosome names separated by newlines).

    optional arguments:
        --params_himut_call                 :   Himut call parameters (default: '"--min_qv 20 --min_mapq 20 --min_bq 20"').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        region_list_file                    :   ${params.region_list_file}
        params_himut_call                   :   ${params_himut_call}
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
workflow VARIANT_CALLING_HIMUT {
    take:
        input_bam_files_ch              // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        region_list_file
        params_himut_call
        output_dir

    main:
        runSamtoolsFaidxFasta(reference_genome_fasta_file)
        fasta_file      = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file  = runSamtoolsFaidxFasta.out.fai_file

        runHimutCall(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            region_list_file,
            params_himut_call,
            output_dir
        )

    emit:
        runHimutCall.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_HIMUT(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.region_list_file,
        params_himut_call,
        params.output_dir
    )
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runPbsv } from '../../modules/pbsv'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.params_pbsv_discover = '--ccs --min-gap-comp-id-perc 97 --min-mapq 20'
params.params_pbsv_call = '--ccs --call-min-reads-per-strand-all-samples 0 --call-min-read-perc-one-sample 10 --call-min-reads-all-samples 3 --call-min-reads-one-sample 3'
params.delete_work_dir = false

if (params.params_pbsv_discover == true) {
    params_pbsv_discover = ''
} else {
    params_pbsv_discover = params.params_pbsv_discover
}

if (params.params_pbsv_call == true) {
    params_pbsv_call = ''
} else {
    params_pbsv_call = params.params_pbsv_call
}

// Step 3. Print inputs and help
log.info """\
         =============================================================================
         Identify structural variants in long-read DNA sequencing BAM files using Pbsv
         =============================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run Pbsv.

    usage: nexus run --nf-workflow long_read_dna_variant_calling_pbsv.nf [required] [optional] [--help]

    required arguments:
        -c                              :   Nextflow .config file.
        -w                              :   Nextflow work directory path.
        --samples_tsv_file              :   TSV file with the following columns:
                                            'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                    :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --params_pbsv_discover          :   Pbsv 'discover' parameters (default: '"--ccs --min-gap-comp-id-perc 97 --min-mapq 20"').
                                            Note that the parameters need to be wrapped in quotes.
        --params_pbsv_call              :   Pbsv 'call' parameters (default:
                                            '"--ccs
                                              --call-min-reads-per-strand-all-samples 0
                                              --call-min-read-perc-one-sample 10
                                              --call-min-reads-all-samples 3
                                              --call-min-reads-one-sample 3"').
                                            Note that the parameters need to be wrapped in quotes.
        --delete_work_dir               :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        reference_genome_fasta_file     :   ${params.reference_genome_fasta_file}
        params_pbsv_discover            :   ${params_pbsv_discover}
        params_pbsv_call                :   ${params_pbsv_call}
        delete_work_dir                 :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow LONG_READ_DNA_VARIANT_CALLING_PBSV {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        params_pbsv_discover
        params_pbsv_call
        output_dir

    main:
        runPbsv(
            input_bam_files_ch,
            reference_genome_fasta_file,
            params_pbsv_discover,
            params_pbsv_call,
            output_dir
        )
}

workflow {
    LONG_READ_DNA_VARIANT_CALLING_PBSV(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params_pbsv_discover,
        params_pbsv_call,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

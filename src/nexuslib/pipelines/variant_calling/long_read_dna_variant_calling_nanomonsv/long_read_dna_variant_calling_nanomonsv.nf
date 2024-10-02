#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runNanomonsv } from '../../modules/nanomonsv'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.params_nanomonsv_parse = ''
params.params_nanomonsv_get = '--min_tumor_variant_read_num 3 --min_tumor_VAF 0.05 --max_control_variant_read_num 0 --max_control_VAF 0.00 --min_indel_size 30 --max_panel_read_num 0 --median_mapQ_thres 20 --qv25'
params.delete_work_dir = false

if (params.params_nanomonsv_parse == true) {
    params_nanomonsv_parse = ''
} else {
    params_nanomonsv_parse = params.params_nanomonsv_parse
}
if (params.params_nanomonsv_get == true) {
    params_nanomonsv_get = ''
} else {
    params_nanomonsv_get = params.params_nanomonsv_get
}

// Step 3. Print inputs and help
log.info """\
         ==========================================================================================
         Identify somatic structural variants in long-read DNA sequencing BAM files using Nanomonsv
         ==========================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run Nanomonsv.

    usage: nexus run --nf-workflow long_read_dna_variant_calling_nanomonsv.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_fasta_fai_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
        --params_nanomonsv_parse            :   Nanomonsv parse parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
        --params_nanomonsv_get              :   Nanomonsv get parameters (default: '"--min_tumor_variant_read_num 3 --min_tumor_VAF 0.05 --max_control_variant_read_num 0 --max_control_VAF 0.00 --min_indel_size 30 --max_panel_read_num 0 --median_mapQ_thres 20 --qv25"').
                                                Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genome_fasta_fai_file     :   ${params.reference_genome_fasta_fai_file}
        params_nanomonsv_parse              :   ${params_nanomonsv_parse}
        params_nanomonsv_get                :   ${params_nanomonsv_get}
        delete_work_dir                     :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.tumor_bam_file}",
        "${row.tumor_bam_bai_file}",
        "${row.normal_bam_file}",
        "${row.normal_bam_bai_file}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow LONG_READ_DNA_VARIANT_CALLING_NANOMONSV {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)]
        reference_genome_fasta_file
        reference_genome_fasta_fai_file
        params_nanomonsv_parse
        params_nanomonsv_get
        output_dir

    main:
        runNanomonsv(
            input_bam_files_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            params_nanomonsv_parse,
            params_nanomonsv_get,
            output_dir
        )
}

workflow {
    LONG_READ_DNA_VARIANT_CALLING_NANOMONSV(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_fai_file,
        params_nanomonsv_parse,
        params_nanomonsv_get,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
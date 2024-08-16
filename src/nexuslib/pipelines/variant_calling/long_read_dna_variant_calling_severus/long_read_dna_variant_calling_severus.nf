#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runSeverus } from '../../modules/severus'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.vntr_bed_file = '/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/severus/Severus-main/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed'
params.params_severus = '--min-support 3 --min-sv-size 30 --min-mapq 20 --output-read-ids --bp-cluster-size 50'
params.delete_work_dir = false

if (params.params_severus == true) {
    params_severus = ''
} else {
    params_severus = params.params_severus
}

// Step 3. Print inputs and help
log.info """\
         ========================================================================================
         Identify somatic structural variants in long-read DNA sequencing BAM files using Severus
         ========================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run Severus.

    usage: nexus run --nf-workflow long_read_dna_variant_calling_severus.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file', 'phased_vcf_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --vntr_bed_file                     :   Tandem repeat regions BED file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/severus/Severus-main/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed).
        --params_severus                    :   Severus parameters (default: '"--min-support 3 --min-sv-size 30 --min-mapq 20 --output-read-ids --bp-cluster-size 50"').
                                                Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        vntr_bed_file                   :   ${params.vntr_bed_file}
        params_severus                  :   ${params_severus}
        delete_work_dir                 :   ${params.delete_work_dir}
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
        "${row.normal_bam_bai_file}",
        "${row.phased_vcf_file}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow LONG_READ_DNA_VARIANT_CALLING_SEVERUS {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        vntr_bed_file
        params_severus
        output_dir

    main:
        runSeverus(
            input_bam_files_ch,
            vntr_bed_file,
            params_severus,
            output_dir
        )
}

workflow {
    LONG_READ_DNA_VARIANT_CALLING_SEVERUS(
        input_bam_files_ch,
        params.vntr_bed_file,
        params_severus,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

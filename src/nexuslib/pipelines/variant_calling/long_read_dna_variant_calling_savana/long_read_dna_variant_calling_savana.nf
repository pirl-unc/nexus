#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runSavanaRun } from '../../modules/savana'
include { runSavanaClassify } from '../../modules/savana'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.custom_params_file = '/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/savana/savana_classification_parameters.json'
params.params_savana_run = '--length 30 --mapq 20 --min_support 3'
params.params_savana_classify = ''
params.delete_work_dir = false

if (params.params_savana_run == true) {
    params_savana_run = ''
} else {
    params_savana_run = params.params_savana_run
}
if (params.params_savana_classify == true) {
    params_savana_classify = ''
} else {
    params_savana_classify = params.params_savana_classify
}

// Step 3. Print inputs and help
log.info """\
         =======================================================================================
         Identify somatic structural variants in long-read DNA sequencing BAM files using Savana
         =======================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run Savana.

    usage: nexus run --nf-workflow long_read_dna_variant_calling_savana.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_fasta_fai_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
        --custom_params_file                :   Savana classify --custom_params file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/savana/savana_classification_parameters.json).
        --params_savana_run                 :   Savana run parameters (default: '"--length 30 --mapq 20 --min_support 3"').
                                                Note that the parameters need to be wrapped in quotes.
        --params_savana_classify            :   Savana classify parameters (default: '""').
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
        custom_params_file                  :   ${params.custom_params_file}
        params_savana_run                   :   ${params_savana_run}
        params_savana_classify              :   ${params_savana_classify}
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
workflow LONG_READ_DNA_VARIANT_CALLING_SAVANA {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)]
        reference_genome_fasta_file
        reference_genome_fasta_fai_file
        custom_params_file
        params_savana_run
        params_savana_classify
        output_dir

    main:
        runSavanaRun(
            input_bam_files_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            params_savana_run,
            output_dir
        )
        runSavanaClassify(
            runSavanaRun.out.f,
            custom_params_file,
            params_savana_classify,
            output_dir
        )
}

workflow {
    LONG_READ_DNA_VARIANT_CALLING_SAVANA(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_fai_file,
        params.custom_params_file,
        params_savana_run,
        params_savana_classify,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
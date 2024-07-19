#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runMantaTumorNormal } from '../../modules/manta'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.params_manta_config = ''
params.params_manta_run = ''
params.delete_work_dir = false

if (params.params_manta_config == true) {
    params_manta_config = ''
} else {
    params_manta_config = params.params_manta_config
}
if (params.params_manta_run == true) {
    params_manta_run = ''
} else {
    params_manta_run = params.params_manta_run
}

// Step 3. Print inputs and help
log.info """\
         =================================================================================
         Identify somatic variants in paired-end read DNA sequencing BAM files using Manta
         =================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run manta (tumor and normal mode).

    usage: nexus run --nf-workflow paired-end_read_dna_variant_calling_manta.nf [required] [optional] [--help]

    required arguments:
        -c                                          :   Nextflow .config file.
        -w                                          :   Nextflow work directory path.
        --samples_tsv_file                          :   TSV file with the following columns:
                                                        'sample_id',
                                                        'tumor_bam_file',
                                                        'tumor_bam_bai_file',
                                                        'normal_bam_file',
                                                        'normal_bam_bai_file'
        --output_dir                                :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file               :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_fasta_fai_file           :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
        --params_manta_config                       :   Manta configManta.py parameters (default: '""').
                                                        Note that the parameters need to be wrapped in quotes.
        --params_manta_run                          :   Manta runWorkflow.py parameters (default: '""').
                                                        Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                           :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                            :   ${params.samples_tsv_file}
        output_dir                                  :   ${params.output_dir}
        reference_genome_fasta_file                 :   ${params.reference_genome_fasta_file}
        reference_genome_fasta_fai_file             :   ${params.reference_genome_fasta_fai_file}
        params_manta_config                         :   ${params_manta_config}
        params_manta_run                            :   ${params_manta_run}
        delete_work_dir                             :   ${params.delete_work_dir}
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
workflow PAIRED_END_READ_DNA_VARIANT_CALLING_MANTA {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(tumor_sample_id), val(normal_sample_id)]
        reference_genome_fasta_file
        reference_genome_fasta_fai_file
        params_manta_config
        params_manta_run
        output_dir

    main:
        runMantaTumorNormal(
            input_bam_files_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            params_manta_config,
            params_manta_run,
            output_dir
        )
}

workflow {
    PAIRED_END_READ_DNA_VARIANT_CALLING_MANTA(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_fai_file,
        params_manta_config,
        params_manta_run,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runDelly2TumorNormal } from '../../modules/delly2'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.exclude_tsv_file = '/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/delly2/human.hg38.excl.tsv'
params.params_delly2call = '--map-qual 20'
params.delete_work_dir = false

if (params.params_delly2call == true) {
    params_delly2call = ''
} else {
    params_delly2call = params.params_delly2call
}

// Step 3. Print inputs and help
log.info """\
         ==================================================================================
         Identify somatic variants in paired-end read DNA sequencing BAM files using Delly2
         ==================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run Delly2 (somatic mode).

    usage: nexus run --nf-workflow paired-end_read_dna_variant_calling_delly2.nf [required] [optional] [--help]

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
        --exclude_tsv_file                          :   Delly2 'call' --exclude TSV file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/delly2/human.hg38.excl.tsv).
        --params_delly2call                         :   Delly2 'call' parameters (default: '"--map-qual 20"').
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
        exclude_tsv_file                            :   ${params.exclude_tsv_file}
        params_delly2call                           :   ${params_delly2call}
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
workflow PAIRED_END_READ_DNA_VARIANT_CALLING_DELLY2 {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(tumor_sample_id), val(normal_sample_id)]
        output_dir
        reference_genome_fasta_file
        reference_genome_fasta_fai_file
        exclude_tsv_file
        params_delly2call

    main:
        runDelly2TumorNormal(
            input_bam_files_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            exclude_tsv_file,
            params_delly2call,
            output_dir
        )
}

workflow {
    PAIRED_END_READ_DNA_VARIANT_CALLING_DELLY2(
        input_bam_files_ch,
        params.output_dir,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_fai_file,
        params.exclude_tsv_file,
        params_delly2call
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

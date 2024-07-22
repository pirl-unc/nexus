#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runStar } from '../../modules/star'
include { copyBamFile } from '../../modules/utils'

// Step 2. Input parameters
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.star_index = '/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/star/hg38_149bp_overhang/'
params.params_star = '--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic'
params.delete_work_dir = false

if (params.params_star == true) {
    params_star = ''
} else {
    params_star = params.params_star
}

// Step 3. Print inputs and help
log.info """\
         ======================================================
         Align paired-end RNA sequencing fastq files using STAR
         ======================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Align paired-end reads to a reference genome index using STAR.

    usage: nexus run --nf-workflow paired-end_read_rna_alignment_star.nf [required] [optional] [--help]

    required arguments:
        -c                              :   Nextflow .config file.
        -w                              :   Nextflow work directory path.
        --samples_tsv_file              :   TSV file with the following columns:
                                            'sample_id', 'fastq_file_1', 'fastq_file_2'.
        --output_dir                    :   Directory to which output files will be copied.

    optional arguments:
        --star_index                    :   Reference genome STAR index (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/star/hg38_100bp_overhang/).
        --params_star                   :   STAR parameters (default: "--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic").
        --delete_work_dir               :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        star_index                      :   ${params.star_index}
        params_star                     :   ${params_star}
        delete_work_dir                 :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.fastq_file_1}",
        "${row.fastq_file_2}") }
    .set { input_fastq_files_ch }

// Step 5. Workflow
workflow PAIRED_END_RNA_READ_ALIGNMENT_STAR {
    take:
        input_fastq_files_ch
        star_index
        params_star
        output_dir
    main:
        runStar(
            input_fastq_files_ch,
            star_index,
            params_star,
            output_dir
        )
        runStar.out.f.set{ copy_bam_file_input_ch }
        copyBamFile(
            copy_bam_file_input_ch,
            output_dir
        )
    emit:
        copyBamFile.out.f
}

workflow {
    PAIRED_END_RNA_READ_ALIGNMENT_STAR(
        input_fastq_files_ch,
        params.star_index,
        params_star,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
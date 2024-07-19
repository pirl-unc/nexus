#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runUltra } from '../../modules/ultra'
include { runSamtoolsCalmd } from '../../modules/samtools'
include { runSamtoolsSort } from '../../modules/samtools'
include { copyBamFile } from '../../modules/utils'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.ultra_index = '/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/ultra/hg38_index/'
params.params_ultra = '--isoseq'
params.delete_work_dir = false

if (params.params_ultra == true) {
    params_ultra = ''
} else {
    params_ultra = params.params_ultra
}

// Step 3. Print inputs and help
log.info """\
         ======================================================
         Align long-read RNA sequencing fastq files using uLTRA
         ======================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Align reads (fastq.gz files) to a reference genome using uLTRA.
        2. Generate MD tags.
        3. Sort MD-tagged bam file.

    usage: nexus run --nf-workflow long_read_rna_alignment_ultra.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'fastq_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
        --ultra_index                       :   uLTRA index path (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/ultra/hg38_index/).
        --params_ultra                      :   uLTRA parameters (default: "--isoseq").
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genome_fasta_fai_file     :   ${params.reference_genome_fasta_fai_file}
        ultra_index                         :   ${params.ultra_index}
        params_ultra                        :   ${params.params_ultra}
        delete_work_dir                     :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.fastq_file}") }
    .set { input_fastq_files_ch }

// Step 5. Workflow
workflow LONG_READ_RNA_ALIGNMENT_ULTRA {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        reference_genome_fasta_file
        reference_genome_fasta_fai_file
        ultra_index
        params_ultra
        output_dir
    main:
        run_ultra_input_ch = input_fastq_files_ch
        runUltra(
            run_ultra_input_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            ultra_index,
            params_ultra,
            output_dir
        )
        runSamtoolsCalmd(
            runUltra.out.f,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file
        )
        runSamtoolsSort(
            runSamtoolsCalmd.out.f
        )
        runSamtoolsSort.out.f.set{ copy_bam_file_input_ch }
        copyBamFile(
            copy_bam_file_input_ch,
            output_dir
        )
    emit:
        copyBamFile.out.f
}

workflow {
    LONG_READ_RNA_ALIGNMENT_ULTRA(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_fai_file,
        params.ultra_index,
        params_ultra,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
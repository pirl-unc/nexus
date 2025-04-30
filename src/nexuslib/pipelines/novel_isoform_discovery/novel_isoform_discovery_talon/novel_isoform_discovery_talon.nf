//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runTalon } from '../../modules/talon'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.reference_genes_gtf_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf'
params.params_talon_init_db = ''
params.params_talon = '--create_novel_spliced_genes'
params.delete_work_dir = false

if (params.params_talon_init_db == true) {
    params_talon_init_db = ''
} else {
    params_talon_init_db = params.params_talon_init_db
}

if (params.params_talon == true) {
    params_talon = ''
} else {
    params_talon = params.params_talon
}

// Step 3. Print inputs and help
log.info """\
         ===================================
         Discover novel isoforms using Talon
         ===================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run Talon.

    usage: nexus run --nf-workflow novel_isoform_discovery_talon.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
        --reference_genes_gtf_file          :   Reference genes GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf).
        --params_talon_init_db              :   talon_initialize_database parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
        --params_talon                      :   talon parameters (default: '"--create_novel_spliced_genes"').
                                                Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genome_fasta_fai_file     :   ${params.reference_genome_fasta_fai_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_talon_init_db                :   ${params_talon_init_db}
        params_talon                        :   ${params_talon}
        delete_work_dir                     :   ${params.delete_work_dir}
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

// Step 4. Workflow
workflow NOVEL_ISOFORM_DISCOVERY_TALON {
    take:
        input_bam_files_ch
        reference_genome_fasta_file
        reference_genome_fasta_fai_file
        reference_genes_gtf_file
        params_talon_init_db
        params_talon
        output_dir
    main:
        runTalon(
            input_bam_files_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            reference_genes_gtf_file,
            params_talon_init_db,
            params_talon,
            output_dir
        )
    emit:
        runTalon.out.f
}

workflow {
    NOVEL_ISOFORM_DISCOVERY_TALON(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_fai_file,
        params.reference_genes_gtf_file,
        params_talon_init_db,
        params_talon,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

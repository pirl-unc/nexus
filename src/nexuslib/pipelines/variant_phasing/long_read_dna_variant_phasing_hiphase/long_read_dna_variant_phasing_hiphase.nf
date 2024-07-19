#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runHiPhaseWith2VcfFiles } from '../../modules/hiphase'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         =======================================================================================
         Phase small and structural variants in long-read DNA sequencing BAM files using HiPhase
         =======================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run HiPhase.

    usage: nexus run --nf-workflow long_read_dna_variant_phasing_hiphase.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id',
                                                'bam_file',
                                                'bam_bai_file',
                                                'small_variants_vcf_gz_file',
                                                'small_variants_vcf_gz_tbi_file',
                                                'structural_variants_vcf_gz_file',
                                                'structural_variants_vcf_gz_tbi_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_fasta_fai_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genome_fasta_fai_file     :   ${params.reference_genome_fasta_fai_file}
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
        "${row.bam_bai_file}",
        "${row.small_variants_vcf_gz_file}",
        "${row.small_variants_vcf_gz_tbi_file}",
        "${row.structural_variants_vcf_gz_file}",
        "${row.structural_variants_vcf_gz_tbi_file}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow LONG_READ_DNA_VARIANT_PHASING_HIPHASE {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file), path(small_variants_vcf_gz_file), path(small_variants_vcf_gz_tbi_file), path(structural_variants_vcf_gz_file), path(structural_variants_vcf_gz_tbi_file)]
        reference_genome_fasta_file
        reference_genome_fasta_fai_file
        output_dir

    main:
        runHiPhaseWith2VcfFiles(
            input_bam_files_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            output_dir
        )
}

workflow {
    LONG_READ_DNA_VARIANT_PHASING_HIPHASE(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_fai_file,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
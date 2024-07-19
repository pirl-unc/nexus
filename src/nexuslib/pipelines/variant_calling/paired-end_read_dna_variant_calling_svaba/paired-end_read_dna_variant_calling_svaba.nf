#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runSvaba } from '../../modules/svaba'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_amb_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.amb'
params.reference_genome_fasta_ann_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.ann'
params.reference_genome_fasta_bwt_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.bwt'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.reference_genome_fasta_pac_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.pac'
params.reference_genome_fasta_sa_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.sa'
params.params_svaba = '--hp --read-tracking'
params.delete_work_dir = false

if (params.params_svaba == true) {
    params_svaba = ''
} else {
    params_svaba = params.params_svaba
}

// Step 3. Print inputs and help
log.info """\
         =================================================================================
         Identify somatic variants in paired-end read DNA sequencing BAM files using Svaba
         =================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run Svaba.

    usage: nexus run --nf-workflow paired-end_read_dna_variant_calling_svaba.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id',
                                                'tumor_bam_file',
                                                'tumor_bam_bai_file',
                                                'normal_bam_file',
                                                'normal_bam_bai_file'
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_fasta_amb_file   :   Reference genome FASTA.AMB file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.amb).
        --reference_genome_fasta_ann_file   :   Reference genome FASTA.ANN file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.ann).
        --reference_genome_fasta_bwt_file   :   Reference genome FASTA.BWT file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.bwt.2bit.64).
        --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
        --reference_genome_fasta_pac_file   :   Reference genome FASTA.PAC file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.pac).
        --reference_genome_fasta_sa_file    :   Reference genome FASTA.PAC file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.sa).
        --params_svaba                      :   Svaba parameters (default: '"--hp --read-tracking"').
                                                Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genome_fasta_amb_file     :   ${params.reference_genome_fasta_amb_file}
        reference_genome_fasta_ann_file     :   ${params.reference_genome_fasta_ann_file}
        reference_genome_fasta_bwt_file     :   ${params.reference_genome_fasta_bwt_file}
        reference_genome_fasta_fai_file     :   ${params.reference_genome_fasta_fai_file}
        reference_genome_fasta_pac_file     :   ${params.reference_genome_fasta_pac_file}
        reference_genome_fasta_sa_file      :   ${params.reference_genome_fasta_sa_file}
        params_svaba                        :   ${params_svaba}
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
workflow PAIRED_END_READ_DNA_VARIANT_CALLING_SVABA {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(tumor_sample_id), val(normal_sample_id)]
        output_dir
        reference_genome_fasta_file
        reference_genome_fasta_amb_file
        reference_genome_fasta_ann_file
        reference_genome_fasta_bwt_file
        reference_genome_fasta_fai_file
        reference_genome_fasta_pac_file
        reference_genome_fasta_sa_file
        params_svaba

    main:
        runSvaba(
            input_bam_files_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_amb_file,
            reference_genome_fasta_ann_file,
            reference_genome_fasta_bwt_file,
            reference_genome_fasta_fai_file,
            reference_genome_fasta_pac_file,
            reference_genome_fasta_sa_file,
            params_svaba,
            output_dir
        )
}

workflow {
    PAIRED_END_READ_DNA_VARIANT_CALLING_SVABA(
        input_bam_files_ch,
        params.output_dir,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_amb_file,
        params.reference_genome_fasta_ann_file,
        params.reference_genome_fasta_bwt_file,
        params.reference_genome_fasta_fai_file,
        params.reference_genome_fasta_pac_file,
        params.reference_genome_fasta_sa_file,
        params_svaba
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

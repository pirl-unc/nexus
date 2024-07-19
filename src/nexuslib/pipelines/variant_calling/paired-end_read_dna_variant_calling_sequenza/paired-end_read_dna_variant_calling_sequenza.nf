#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runSequenzaUtilsBam2Seqz } from '../../modules/sequenza'
include { mergeSequenzaSeqzFiles } from '../../modules/sequenza'
include { runSequenzaUtilsSeqzBinning } from '../../modules/sequenza'
include { runSequenza } from '../../modules/sequenza'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''

// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_wig_file = '/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/sequenza/hg38.gc50.wig.gz'
params.assembly = 'hg38'
params.chromosomes = 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 ch22 chrX chrY'
params.params_sequenza_bam2seqz = '-N 20 --qformat sanger'
params.params_sequenza_seqzbinning = '--window 50'
params.delete_work_dir = false

if (params.params_sequenza_bam2seqz == true) {
    params_sequenza_bam2seqz = ''
} else {
    params_sequenza_bam2seqz = params.params_sequenza_bam2seqz
}

if (params.params_sequenza_seqzbinning == true) {
    params_sequenza_seqzbinning = ''
} else {
    params_sequenza_seqzbinning = params.params_sequenza_seqzbinning
}

// Step 3. Print inputs and help
log.info """\
         ================================================================================================
         Identify somatic copy number variants in paired-end read DNA sequencing BAM files using Sequenza
         ================================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1.  Run sequenza-utils bam2seqz.
        2.  Merge seqz files into one seqz file.
        3.  Run sequenza-utils seqz_binning.
        4.  Run sequenza in R.

    usage: nexus run --nf-workflow paired-end_read_dna_variant_calling_sequenza.nf [required] [optional] [--help]

    required arguments:
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_wig_file         :   Reference genome WIG file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/sequenza/hg38.gc50.wig.gz).
        --assembly                          :   Assembly ('hg19' or 'hg38').
        --chromosomes                       :   List of chromosomes (e.g. 'chr1 chr2 chr3').
        --params_sequenza_bam2seqz          :   Sequenza 'bam2seqz' parameters (default: '"-N 20 --qformat sanger"').
                                                Note that the parameters need to be wrapped in quotes.
        --params_sequenza_seqzbinning       :   Sequenza 'seqz_binning' parameters (default: '"--window 50"').
                                                Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genome_wig_file           :   ${params.reference_genome_wig_file}
        assembly                            :   ${params.assembly}
        chromosomes                         :   ${params.chromosomes}
        params_sequenza_bam2seqz            :   ${params_sequenza_bam2seqz}
        params_sequenza_seqzbinning         :   ${params_sequenza_seqzbinning}
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

// Step 5. Workflows
workflow PAIRED_END_READ_DNA_VARIANT_CALLING_SEQUENZA {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(tumor_sample_id), val(normal_sample_id)]
        output_dir
        reference_genome_fasta_file
        reference_genome_wig_file
        assembly
        chromosomes
        params_sequenza_bam2seqz
        params_sequenza_seqzbinning

    main:
        runSequenzaUtilsBam2Seqz(
            input_bam_files_ch,
            reference_genome_fasta_file,
            reference_genome_wig_file,
            chromosomes,
            params_sequenza_bam2seqz,
            output_dir
        )
        mergeSequenzaSeqzFiles(
            runSequenzaUtilsBam2Seqz.out.f,
            output_dir
        )
        runSequenzaUtilsSeqzBinning(
            mergeSequenzaSeqzFiles.out.f,
            params_sequenza_seqzbinning,
            output_dir
        )
        chromosomes = chromosomes.replace(" ",",")
        runSequenza(
            runSequenzaUtilsSeqzBinning.out.f,
            chromosomes,
            assembly,
            output_dir
        )
}

workflow {
    PAIRED_END_READ_DNA_VARIANT_CALLING_SEQUENZA(
        input_bam_files_ch,
        params.output_dir,
        params.reference_genome_fasta_file,
        params.reference_genome_wig_file,
        params.assembly,
        params.chromosomes,
        params_sequenza_bam2seqz,
        params_sequenza_seqzbinning
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }                from '../../../tools/samtools'
include { runSequenzaUtilsIndex }           from '../../../tools/sequenza'
include { runSequenzaUtilsBam2Seqz }        from '../../../tools/sequenza'
include { runSequenzaUtilsSeqzBinning }     from '../../../tools/sequenza'
include { mergeSequenzaSeqzFiles }          from '../../../tools/sequenza'
include { runSequenza }                     from '../../../tools/sequenza'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help = ''

// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.reference_genome_fasta_file      = ''
params.assembly                         = ''

// Optional arguments
params.sequenzautils_gcwiggle           = '-w 50'
params.sequenzautils_bam2seqz           = '-N 20 --qformat sanger'
params.sequenzautils_seqzbinning        = '--window 50'
params.chromosomes                      = 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 ch22 chrX chrY'

if (params.sequenzautils_gcwiggle == true) {
    params_sequenzautils_gcwiggle = ''
} else {
    params_sequenzautils_gcwiggle = params.sequenzautils_gcwiggle
}

if (params.sequenzautils_bam2seqz == true) {
    params_sequenzautils_bam2seqz = ''
} else {
    params_sequenzautils_bam2seqz = params.sequenzautils_bam2seqz
}

if (params.sequenzautils_seqzbinning == true) {
    params_sequenzautils_seqzbinning = ''
} else {
    params_sequenzautils_seqzbinning = params.sequenzautils_seqzbinning
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
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
                                                'sample_id',
                                                'tumor_bam_file',
                                                'tumor_bam_bai_file',
                                                'normal_bam_file',
                                                'normal_bam_bai_file'
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.
        --assembly                          :   Assembly ('hg19' or 'hg38').

    optional arguments:
        --chromosomes                       :   List of chromosomes (default: '"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 ch22 chrX chrY"').
                                                Note that the parameters need to be wrapped in quotes.
        --params_sequenzautils_gcwiggle     :   sequenza-utils 'gc_wiggle' parameters (default: '"-w 50"').
                                                Note that the parameters need to be wrapped in quotes.
        --params_sequenzautils_bam2seqz     :   sequenza-utils 'bam2seqz' parameters (default: '"-N 20 --qformat sanger"').
                                                Note that the parameters need to be wrapped in quotes.
        --params_sequenzautils_seqzbinning  :   sequenza-utils 'seqz_binning' parameters (default: '"--window 50"').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        assembly                            :   ${params.assembly}
        chromosomes                         :   ${params.chromosomes}
        params_sequenzautils_gcwiggle       :   ${params_sequenzautils_gcwiggle}
        params_sequenzautils_bam2seqz       :   ${params_sequenzautils_bam2seqz}
        params_sequenzautils_seqzbinning    :   ${params_sequenzautils_seqzbinning}
    """.stripIndent()
}

// ------------------------------------------------------------
// Step 4. Set channels
// ------------------------------------------------------------
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

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_SEQUENZA {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)]
        reference_genome_fasta_file
        assembly
        chromosomes
        params_sequenzautils_gcwiggle
        params_sequenzautils_bam2seqz
        params_sequenzautils_seqzbinning
        output_dir

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file                  = runSamtoolsFaidx.out.fasta
        fasta_fai_file              = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file              = runSamtoolsFaidx.out.gzi_file

        runSequenzaUtilsIndex(
            fasta_file,
            params_sequenzautils_gcwiggle
        )
        gc_wiggle_file = runSequenzaUtilsIndex.out.wig_file

        runSequenzaUtilsBam2Seqz(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            fasta_gzi_file,
            gc_wiggle_file,
            chromosomes,
            params_sequenzautils_bam2seqz,
            output_dir
        )

        mergeSequenzaSeqzFiles(
            runSequenzaUtilsBam2Seqz.out.f,
            output_dir
        )

        runSequenzaUtilsSeqzBinning(
            mergeSequenzaSeqzFiles.out.f,
            params_sequenzautils_seqzbinning,
            output_dir
        )

        chromosomes = chromosomes.replace(" ",",")

        runSequenza(
            runSequenzaUtilsSeqzBinning.out.f,
            chromosomes,
            assembly,
            output_dir
        )

    emit:
        runSequenza.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_SEQUENZA(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.assembly,
        params.chromosomes,
        params_sequenzautils_gcwiggle,
        params_sequenzautils_bam2seqz,
        params_sequenzautils_seqzbinning,
        params.output_dir
    )
}

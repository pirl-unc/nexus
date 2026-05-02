#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }                  from '../../../tools/samtools'
include { runLongphasePhase }                      from '../../../tools/longphase'
include { runLongphaseHaplotag }                   from '../../../tools/longphase'
include { decompressFile as decompressFasta }      from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''

// Optional arguments
// `longphase phase` v2.0.1 requires exactly one of '--ont' or '--pb'
// (PacBio HiFi/CCS data uses '--pb'). `longphase haplotag` has no platform
// flag — leave params_longphase_haplotag empty unless you want to pass
// non-platform options.
params.params_longphase_phase           = '--ont'
params.params_longphase_haplotag        = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow HAPLOTAGGING_LONGPHASE {
    take:
        input_bam_vcf_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file), path(small_variants_vcf_file)]
        reference_genome_fasta_file
        params_longphase_phase
        params_longphase_haplotag
        output_dir
        subdir                             // sub-folder appended after ${output_dir}/${sample_id}/; pass '' for no nesting.

    main:
        // Step 1. Decompress and index reference genome FASTA file
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidxFasta(decompressFasta.out.f)
        fasta_file          = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file      = runSamtoolsFaidxFasta.out.fai_file

        // Step 2. Run LongPhase 'phase' to produce a phased VCF
        runLongphasePhase(
            input_bam_vcf_files_ch,
            fasta_file,
            fasta_fai_file,
            params_longphase_phase,
            output_dir,
            subdir
        )

        // Step 3. Join phased VCF (and its index) with the original BAM channel by sample_id,
        //         producing the input expected by runLongphaseHaplotag:
        //         [sample_id, bam, bai, phased_vcf, phased_vcf_tbi]
        bam_only_ch = input_bam_vcf_files_ch
            .map { sample_id, bam_file, bam_bai_file, vcf_file -> tuple(sample_id, bam_file, bam_bai_file) }

        haplotag_input_ch = bam_only_ch
            .join(runLongphasePhase.out.f, by: 0)
            .map { sample_id, bam_file, bam_bai_file, phased_vcf, phased_vcf_tbi ->
                tuple(sample_id, bam_file, bam_bai_file, phased_vcf, phased_vcf_tbi) }

        // Step 4. Run LongPhase 'haplotag' to produce a haplotagged BAM
        runLongphaseHaplotag(
            haplotag_input_ch,
            fasta_file,
            fasta_fai_file,
            params_longphase_haplotag,
            output_dir,
            subdir
        )

    emit:
        phased_vcf      = runLongphasePhase.out.f         // [sample_id, vcf.gz, vcf.gz.tbi]
        haplotagged_bam = runLongphaseHaplotag.out.f      // [sample_id, bam, bam.bai]
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===============================================================================================
             Phase small variants and haplotag long-read DNA sequencing BAM files using LongPhase
             ===============================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run LongPhase 'phase' command to produce a phased VCF.
            2. Run LongPhase 'haplotag' command to produce a haplotagged BAM (using the phased VCF from step 1).

        usage: nexus run --nf-workflow haplotagging_longphase.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id',
                                                    'bam_file',
                                                    'bam_bai_file',
                                                    'small_variants_vcf_file'
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.

        optional arguments:
            --params_longphase_phase            :   LongPhase 'phase' parameters (default: '"--ont"').
                                                    Set to '--pb' (PacBio HiFi/CCS) or '--ont' (Oxford
                                                    Nanopore) to match the sequencing platform. Note
                                                    that the parameters need to be wrapped in quotes.
            --params_longphase_haplotag         :   LongPhase 'haplo    tag' parameters (default: '""').
                                                    `longphase haplotag` has no platform flag — leave
                                                    empty unless passing other options. Note that the
                                                    parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_longphase_phase    = (params.params_longphase_phase    == true) ? '' : params.params_longphase_phase
    def params_longphase_haplotag = (params.params_longphase_haplotag == true) ? '' : params.params_longphase_haplotag

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_longphase_phase              :   ${params_longphase_phase}
        params_longphase_haplotag           :   ${params_longphase_haplotag}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}",
            "${row.small_variants_vcf_file}") }
        .set { input_bam_vcf_files_ch }

    HAPLOTAGGING_LONGPHASE(
        input_bam_vcf_files_ch,
        params.reference_genome_fasta_file,
        params_longphase_phase,
        params_longphase_haplotag,
        params.output_dir,
        ''     // standalone: don't nest outputs under a subdir
    )
}

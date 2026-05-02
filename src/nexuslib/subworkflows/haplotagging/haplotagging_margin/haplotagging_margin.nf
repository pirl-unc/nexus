#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }                  from '../../../tools/samtools'
include { runMarginPhase }                         from '../../../tools/margin'
include { decompressFile as decompressFasta }      from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''

// Margin requires an HMM parameters JSON file. Required — must be supplied
// by the user. The margin docker image ships defaults under /opt/margin/params/
// (e.g. phase/allParams.phase_vcf.pb-hifi.json) but those paths are inside
// the container; the user must point to a host-accessible file. Test
// fixtures live at test/data/indices/margin/.
//
// NOTE: Margin v2.3.1 dropped the standalone `haplotag` subcommand —
// `margin phase --produceFinalHaplotaggedBam` produces the phased VCF AND
// the haplotagged BAM in a single call, so only the phase params JSON is
// needed here.
params.margin_phase_params_json_file    = ''

// Optional argument string appended verbatim to each margin invocation.
params.params_margin_phase              = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow HAPLOTAGGING_MARGIN {
    take:
        input_bam_vcf_files_ch              // channel: [val(sample_id), path(bam_file), path(bam_bai_file), path(small_variants_vcf_file)]
        reference_genome_fasta_file
        margin_phase_params_json_file       // val: path string to phase params JSON
        params_margin_phase
        output_dir
        subdir                              // sub-folder appended after ${output_dir}/${sample_id}/; pass '' for no nesting.

    main:
        // Step 1. Decompress and index reference genome FASTA file
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidxFasta(decompressFasta.out.f)
        fasta_file          = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file      = runSamtoolsFaidxFasta.out.fai_file

        // Step 2. Run `margin phase --produceFinalHaplotaggedBam` to produce
        //         BOTH the phased VCF and the haplotagged BAM in one call.
        runMarginPhase(
            input_bam_vcf_files_ch,
            fasta_file,
            fasta_fai_file,
            margin_phase_params_json_file,
            params_margin_phase,
            output_dir,
            subdir
        )

    emit:
        phased_vcf      = runMarginPhase.out.phased_vcf        // [sample_id, vcf.gz, vcf.gz.tbi]
        haplotagged_bam = runMarginPhase.out.haplotagged_bam   // [sample_id, bam, bam.bai]
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===============================================================================================
             Phase small variants and haplotag long-read DNA sequencing BAM files using Margin
             ===============================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run `margin phase --produceFinalHaplotaggedBam` to produce a phased VCF
               and a haplotagged BAM in one call (Margin v2.3.1 has no separate
               `haplotag` subcommand).

        usage: nexus run --nf-workflow haplotagging_margin.nf [required] [optional] [--help]

        required arguments:
            -c                                       :   Nextflow .config file.
            -w                                       :   Nextflow work directory path.
            --samples_tsv_file                       :   TSV file with the following columns:
                                                         'sample_id',
                                                         'bam_file',
                                                         'bam_bai_file',
                                                         'small_variants_vcf_file'
            --output_dir                             :   Directory to which output files will be copied.
            --reference_genome_fasta_file            :   Reference genome FASTA file.
            --margin_phase_params_json_file          :   Margin 'phase' HMM parameters JSON file. Pick a file
                                                         matching your sequencing platform — sample fixtures
                                                         ship at test/data/indices/margin/phase/, and the
                                                         margin docker image carries the same set under
                                                         /opt/margin/params/.

        optional arguments:
            --params_margin_phase                    :   Margin 'phase' extra arguments appended verbatim
                                                         (default: '""'). Wrap in quotes.
        """.stripIndent()
        exit 0
    }

    def params_margin_phase = (params.params_margin_phase == true) ? '' : params.params_margin_phase

    log.info"""\
        samples_tsv_file                       :   ${params.samples_tsv_file}
        output_dir                             :   ${params.output_dir}
        reference_genome_fasta_file            :   ${params.reference_genome_fasta_file}
        margin_phase_params_json_file          :   ${params.margin_phase_params_json_file}
        params_margin_phase                    :   ${params_margin_phase}
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

    HAPLOTAGGING_MARGIN(
        input_bam_vcf_files_ch,
        params.reference_genome_fasta_file,
        params.margin_phase_params_json_file,
        params_margin_phase,
        params.output_dir,
        ''     // standalone: don't nest outputs under a subdir
    )
}

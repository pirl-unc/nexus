#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 *
 * Phase variants with HapCUT2 and tag reads with WhatsHap haplotag.
 * HapCUT2 produces only a phased VCF; this subworkflow chains in WhatsHap's
 * haplotag step so users get both the phased VCF and a haplotagged BAM in
 * one call (parallel to the merged HAPLOTAGGING_WHATSHAP design).
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }                  from '../../../tools/samtools'
include { runHapCUT2 }                             from '../../../tools/hapcut2'
include { runWhatshapHaplotag }                    from '../../../tools/whatshap'
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
params.read_technology                  = 'pacbio'   // 'pacbio' | 'ont' | 'illumina'
params.params_extracthairs              = ''
params.params_hapcut2                   = ''
params.params_whatshap_haplotag         = '--ignore-read-groups --skip-missing-contigs --output-threads 4'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow HAPLOTAGGING_HAPCUT2_WHATSHAP {
    take:
        input_bam_vcf_files_ch          // channel: [val(sample_id), path(bam_file), path(bam_bai_file), path(small_variants_vcf_file)]
        reference_genome_fasta_file
        read_technology
        params_extracthairs
        params_hapcut2
        params_whatshap_haplotag
        output_dir
        subdir                          // sub-folder appended after ${output_dir}/${sample_id}/; pass '' for no nesting.

    main:
        // Step 1. Decompress and index reference genome FASTA file
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidxFasta(decompressFasta.out.f)
        fasta_file          = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file      = runSamtoolsFaidxFasta.out.fai_file

        // Step 2. Run HapCUT2 to produce a phased VCF
        runHapCUT2(
            input_bam_vcf_files_ch,
            fasta_file,
            fasta_fai_file,
            read_technology,
            params_extracthairs,
            params_hapcut2,
            output_dir,
            subdir
        )

        // Step 3. Build haplotag input by joining the BAM channel with HapCUT2's
        //         phased VCF + index. runHapCUT2 emits:
        //         (sid, haplotypes.txt, phased.vcf.gz, phased.vcf.gz.tbi)
        //         runWhatshapHaplotag expects:
        //         (sid, bam, bai, phased_vcf, phased_vcf_tbi)
        bam_only_ch = input_bam_vcf_files_ch
            .map { sample_id, bam_file, bam_bai_file, vcf_file ->
                tuple(sample_id, bam_file, bam_bai_file) }

        haplotag_input_ch = bam_only_ch
            .join(runHapCUT2.out.f, by: 0)
            .map { sample_id, bam_file, bam_bai_file, hapblocks, phased_vcf, phased_vcf_tbi ->
                tuple(sample_id, bam_file, bam_bai_file, phased_vcf, phased_vcf_tbi) }

        // Step 4. Run WhatsHap haplotag to produce a haplotagged BAM
        runWhatshapHaplotag(
            haplotag_input_ch,
            fasta_file,
            fasta_fai_file,
            params_whatshap_haplotag,
            output_dir,
            subdir
        )

    emit:
        phased_vcf      = runHapCUT2.out.f               // [sample_id, haplotypes.txt, vcf.gz, vcf.gz.tbi]
        haplotagged_bam = runWhatshapHaplotag.out.f      // [sample_id, bam, bam.bai]
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===============================================================================================
             Phase small variants with HapCUT2 and haplotag long-read DNA sequencing BAM files using WhatsHap
             ===============================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run extractHAIRS + HAPCUT2 to produce a phased VCF.
            2. Run Whatshap 'haplotag' on the phased VCF to produce a haplotagged BAM.

        usage: nexus run --nf-workflow haplotagging_hapcut2-whatshap.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id',
                                                    'bam_file',
                                                    'bam_bai_file',
                                                    'small_variants_vcf_file'
                                                    (VCF may be plain '.vcf' or '.vcf.gz')
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.

        optional arguments:
            --read_technology                   :   Read technology for extractHAIRS
                                                    (choices: pacbio, ont, illumina; default: 'pacbio').
            --params_extracthairs               :   Extra extractHAIRS parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
            --params_hapcut2                    :   Extra HAPCUT2 parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
            --params_whatshap_haplotag          :   Whatshap 'haplotag' parameters (default:
                                                    '"--ignore-read-groups --skip-missing-contigs --output-threads 4"').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_extracthairs       = (params.params_extracthairs       == true) ? '' : params.params_extracthairs
    def params_hapcut2            = (params.params_hapcut2            == true) ? '' : params.params_hapcut2
    def params_whatshap_haplotag  = (params.params_whatshap_haplotag  == true) ? '' : params.params_whatshap_haplotag

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        read_technology                     :   ${params.read_technology}
        params_extracthairs                 :   ${params_extracthairs}
        params_hapcut2                      :   ${params_hapcut2}
        params_whatshap_haplotag            :   ${params_whatshap_haplotag}
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

    HAPLOTAGGING_HAPCUT2_WHATSHAP(
        input_bam_vcf_files_ch,
        params.reference_genome_fasta_file,
        params.read_technology,
        params_extracthairs,
        params_hapcut2,
        params_whatshap_haplotag,
        params.output_dir,
        ''     // standalone: don't nest outputs under a subdir
    )
}

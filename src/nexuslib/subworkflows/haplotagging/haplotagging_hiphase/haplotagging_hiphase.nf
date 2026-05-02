#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }                                  from '../../../tools/samtools'
include { runHiPhaseWith2VcfFiles }                           from '../../../tools/hiphase'
include { decompressFile as decompressFasta }                 from '../../../tools/utils'
include { bgzipAndIndexVcfFile as bgzipAndIndexSmallVcf }     from '../../../tools/utils'
include { bgzipAndIndexVcfFile as bgzipAndIndexStructVcf }    from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow HAPLOTAGGING_HIPHASE {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        small_variants_vcf_ch          // channel: [val(sample_id), path(small_variants_vcf_file)]    — may be plain .vcf or .vcf.gz
        structural_variants_vcf_ch     // channel: [val(sample_id), path(structural_variants_vcf_file)] — may be plain .vcf or .vcf.gz
        reference_genome_fasta_file
        output_dir
        subdir                         // sub-folder appended after ${output_dir}/${sample_id}/; pass '' for no nesting.

    main:
        // Step 1. Decompress and re-bgzip+index reference genome FASTA file
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidx(decompressFasta.out.f)
        fasta_file          = runSamtoolsFaidx.out.fasta
        fasta_fai_file      = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file      = runSamtoolsFaidx.out.gzi_file

        // Step 2. Ensure both VCFs are bgzipped and tabix-indexed
        bgzipAndIndexSmallVcf(small_variants_vcf_ch)
        bgzipAndIndexStructVcf(structural_variants_vcf_ch)

        // Step 3. Join indexed VCFs with the BAM channel by sample_id, producing the
        //         input expected by runHiPhaseWith2VcfFiles:
        //         [sample_id, bam, bai, small_vcf.gz, small_vcf.gz.tbi, struct_vcf.gz, struct_vcf.gz.tbi]
        hiphase_input_ch = input_bam_files_ch
            .join(bgzipAndIndexSmallVcf.out.f, by: 0)
            .join(bgzipAndIndexStructVcf.out.f, by: 0)

        // Step 4. Run HiPhase
        runHiPhaseWith2VcfFiles(
            hiphase_input_ch,
            fasta_file,
            fasta_fai_file,
            output_dir,
            subdir
        )

    emit:
        runHiPhaseWith2VcfFiles.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =======================================================================================
             Phase small and structural variants in long-read DNA sequencing BAM files using HiPhase
             =======================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. bgzip and tabix-index the small and structural variant VCFs (skipped if already .gz).
            2. Run HiPhase.

        usage: nexus run --nf-workflow haplotagging_hiphase.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id',
                                                    'bam_file',
                                                    'bam_bai_file',
                                                    'small_variants_vcf_file',
                                                    'structural_variants_vcf_file'.
                                                    VCF files may be either plain '.vcf' or '.vcf.gz';
                                                    they will be bgzipped and indexed automatically.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
        """.stripIndent()
        exit 0
    }

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
    """.stripIndent()

    // BAM channel
    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    // Small variants VCF channel
    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.small_variants_vcf_file}") }
        .set { small_variants_vcf_ch }

    // Structural variants VCF channel
    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.structural_variants_vcf_file}") }
        .set { structural_variants_vcf_ch }

    HAPLOTAGGING_HIPHASE(
        input_bam_files_ch,
        small_variants_vcf_ch,
        structural_variants_vcf_ch,
        params.reference_genome_fasta_file,
        params.output_dir,
        ''     // standalone: don't nest outputs under a subdir
    )
}

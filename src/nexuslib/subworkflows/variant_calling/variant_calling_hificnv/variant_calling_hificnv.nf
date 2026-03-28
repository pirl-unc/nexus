#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runHiFiCNV }                             from '../../../tools/hificnv'
include { decompressFile as decompressFasta }      from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.exclude_bed_file                 = ''

// Optional arguments
params.params_hificnv                   = ''

if (params.params_hificnv == true) {
    params_hificnv = ''
} else {
    params_hificnv = params.params_hificnv
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ====================================================================================
         Identify copy number alterations in long-read DNA sequencing BAM files using HiFiCNV
         ====================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run hificnv.

    usage: nexus run --nf-workflow variant_calling_hificnv.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'bam_file', 'bam_bai_file', 'vcf_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.
        --exclude_bed_file                  :   Exclude BED file.

    optional arguments:
        --params_hificnv                    :   hificnv parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        exclude_bed_file                    :   ${params.exclude_bed_file}
        params_hificnv                      :   ${params_hificnv}
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
        "${row.bam_file}",
        "${row.bam_bai_file}",
        "${row.vcf_file}") }
    .set { input_bam_vcf_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_HIFICNV {
    take:
        input_bam_vcf_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file), path(vcf_file)]
        reference_genome_fasta_file
        exclude_bed_file
        params_hificnv
        output_dir

    main:
        decompressFasta(reference_genome_fasta_file)

        runHiFiCNV(
            input_bam_vcf_files_ch,
            decompressFasta.out.f,
            exclude_bed_file,
            params_hificnv,
            output_dir
        )

    emit:
        runHiFiCNV.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_HIFICNV(
        input_bam_vcf_files_ch,
        params.reference_genome_fasta_file,
        params.exclude_bed_file,
        params_hificnv,
        params.output_dir
    )
}

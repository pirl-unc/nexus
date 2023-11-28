#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runDeepVariantSingularity } from '../../modules/deepvariant'
include { runDeepVariantDocker } from '../../modules/deepvariant'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.reference_genome_fasta_file = ''
params.containerization = 'singularity' // or docker
params.singularity = ''
params.deepvariant_bin_path = ''
params.deepvariant_bin_version = ''
params.deepvariant_lib_path = ''
params.deepvariant_output_path = ''
params.deepvariant_model_type = ''
// Optional arguments
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         ==================================================================================
         Identify Small Variants (SNVs and INDELs) using Long-read DNA Sequencing BAM Files
         ==================================================================================
         PURPOSE:
            This workflow is intended for long-read DNA sequencing BAM files.

         WORKFLOW:
            1. Run DeepVariant.
         """.stripIndent()

if (params.help) {
    log.info"""\
        usage: nextflow run long_read_dna_small_variants.nf [required] [optional] [--help]

        required arguments:
            --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'bam_file', 'bam_bai_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --singularity                       :   Singularity path.
            --deepvariant_bin_path              :   DeepVariant bin path.
            --deepvariant_bin_version           :   DeepVariant bin version.
            --deepvariant_lib_path              :   DeepVariant lib path.
            --deepvariant_output_path           :   DeepVariant output path.
            --deepvariant_model_type            :   DeepVariant --model_type parameter value.

        optional arguments:
            --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
            samples_tsv_file                    :   ${params.samples_tsv_file}
            output_dir                          :   ${params.output_dir}
            reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
            containerization                    :   ${params.containerization}
            deepvariant_bin_path                :   ${params.deepvariant_bin_path}
            deepvariant_bin_version             :   ${params.deepvariant_bin_version}
            deepvariant_lib_path                :   ${params.deepvariant_lib_path}
            deepvariant_output_path             ::  ${params.deepvariant_output_path}
            deepvariant_model_type              :   ${params.deepvariant_model_type}
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

// Step 5. Workflow
workflow PACBIO_DNA_SMALL_VARIANTS {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        containerization
        deepvariant_bin_version
        deepvariant_bin_path
        deepvariant_lib_path
        deepvariant_output_path
        deepvariant_model_type
        output_dir

    main:
        input_bam_files_ch
            .map{ [it[0], it[1], it[2]] }
            .set{ run_deepvariant_input_ch }

        // DeepVariant
        if (containerization == "singularity") {
            runDeepVariantSingularity(
                run_deepvariant_input_ch,
                reference_genome_fasta_file,
                containerization,
                deepvariant_model_type,
                deepvariant_bin_version,
                deepvariant_bin_path,
                deepvariant_lib_path,
                output_dir
            )
        }
        if (containerization == "docker") {
            runDeepVariantDocker(
                run_deepvariant_input_ch,
                reference_genome_fasta_file,
                containerization,
                deepvariant_model_type,
                deepvariant_bin_version,
                deepvariant_bin_path,
                deepvariant_lib_path,
                deepvariant_output_path,
                output_dir
            )
        }
}

workflow {
    PACBIO_DNA_SMALL_VARIANTS(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.containerization,
        params.deepvariant_bin_version,
        params.deepvariant_bin_path,
        params.deepvariant_lib_path,
        params.deepvariant_output_path,
        params.deepvariant_model_type,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

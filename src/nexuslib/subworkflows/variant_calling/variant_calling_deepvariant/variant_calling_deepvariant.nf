#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }              from '../../../tools/samtools'
include { runDeepVariantSingularity }          from '../../../tools/deepvariant'
include { runDeepVariantDocker }               from '../../../tools/deepvariant'
include { decompressFile as decompressFasta }  from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.deepvariant_input_path           = ''
params.deepvariant_output_path          = ''

// Optional arguments
params.deepvariant_containerization     = 'singularity' // or docker
params.deepvariant_model_type           = 'PACBIO'
params.deepvariant_bin_path             = '/opt/deepvariant/bin/run_deepvariant'
params.deepvariant_bin_version          = '1.9.0'

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         =================================================================================================
         Identify small variants (SNVs and INDELs) in long-read DNA sequencing BAM files using DeepVariant
         =================================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run DeepVariant.

    usage: nexus run --nf-workflow variant_calling_deepvariant.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.
        --deepvariant_input_path            :   DeepVariant input path.
        --deepvariant_output_path           :   DeepVariant output path.

    optional arguments:
        --deepvariant_containerization      :   Containerization ('singularity' or 'docker'; default: 'singularity').
        --deepvariant_model_type            :   DeepVariant --model_type parameter value (default: 'PACBIO').
        --deepvariant_bin_path              :   DeepVariant bin path (default: '/opt/deepvariant/bin/run_deepvariant').
        --deepvariant_bin_version           :   DeepVariant bin version (default: '1.9.0').
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        deepvariant_containerization        :   ${params.deepvariant_containerization}
        deepvariant_bin_path                :   ${params.deepvariant_bin_path}
        deepvariant_bin_version             :   ${params.deepvariant_bin_version}
        deepvariant_input_path              :   ${params.deepvariant_input_path}
        deepvariant_output_path             :   ${params.deepvariant_output_path}
        deepvariant_model_type              :   ${params.deepvariant_model_type}
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
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_DEEPVARIANT {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        deepvariant_containerization
        deepvariant_bin_version
        deepvariant_bin_path
        deepvariant_input_path
        deepvariant_output_path
        deepvariant_model_type
        output_dir

    main:
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidxFasta(decompressFasta.out.f)
        fasta_file      = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file  = runSamtoolsFaidxFasta.out.fai_file

        input_bam_files_ch
            .map{ [it[0], it[1], it[2]] }
            .set{ run_deepvariant_input_ch }

        // DeepVariant
        if (deepvariant_containerization == "singularity") {
            runDeepVariantSingularity(
                run_deepvariant_input_ch,
                fasta_file,
                fasta_fai_file,
                deepvariant_containerization,
                deepvariant_model_type,
                deepvariant_bin_version,
                deepvariant_bin_path,
                deepvariant_input_path,
                deepvariant_output_path,
                output_dir
            )
        }
        if (deepvariant_containerization == "docker") {
            runDeepVariantDocker(
                run_deepvariant_input_ch,
                fasta_file,
                fasta_fai_file,
                deepvariant_model_type,
                deepvariant_bin_version,
                deepvariant_bin_path,
                deepvariant_input_path,
                deepvariant_output_path,
                output_dir
            )
        }

        // Collect output from whichever containerization was used
        deepvariant_out_ch = Channel.empty()
        if (deepvariant_containerization == "singularity") {
            deepvariant_out_ch = runDeepVariantSingularity.out.f
        }
        if (deepvariant_containerization == "docker") {
            deepvariant_out_ch = runDeepVariantDocker.out.f
        }

    emit:
        deepvariant_out_ch
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_DEEPVARIANT(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.deepvariant_containerization,
        params.deepvariant_bin_version,
        params.deepvariant_bin_path,
        params.deepvariant_input_path,
        params.deepvariant_output_path,
        params.deepvariant_model_type,
        params.output_dir
    )
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }              from '../../../tools/samtools'
include { runDeepSomaticSingularity }          from '../../../tools/deepsomatic'
include { runDeepSomaticDocker }               from '../../../tools/deepsomatic'
include { decompressFile as decompressFasta }  from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.deepsomatic_input_path           = ''
params.deepsomatic_output_path          = ''

// Optional arguments
params.deepsomatic_containerization     = 'singularity' // or docker
params.deepsomatic_model_type           = 'PACBIO'
params.deepsomatic_bin_path             = 'run_deepsomatic'
params.deepsomatic_bin_version          = '1.9.0'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_DEEPSOMATIC {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)]
        reference_genome_fasta_file
        deepsomatic_containerization
        deepsomatic_bin_version
        deepsomatic_bin_path
        deepsomatic_input_path
        deepsomatic_output_path
        deepsomatic_model_type
        output_dir

    main:
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidxFasta(decompressFasta.out.f)
        fasta_file      = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file  = runSamtoolsFaidxFasta.out.fai_file

        if (deepsomatic_containerization == "singularity") {
            runDeepSomaticSingularity(
                input_bam_files_ch,
                fasta_file,
                fasta_fai_file,
                deepsomatic_containerization,
                deepsomatic_model_type,
                deepsomatic_bin_version,
                deepsomatic_bin_path,
                deepsomatic_input_path,
                deepsomatic_output_path,
                output_dir
            )
        }
        if (deepsomatic_containerization == "docker") {
            runDeepSomaticDocker(
                input_bam_files_ch,
                fasta_file,
                fasta_fai_file,
                deepsomatic_model_type,
                deepsomatic_bin_version,
                deepsomatic_bin_path,
                deepsomatic_input_path,
                deepsomatic_output_path,
                output_dir
            )
        }

        // Collect output from whichever containerization was used
        deepsomatic_out_ch = Channel.empty()
        if (deepsomatic_containerization == "singularity") {
            deepsomatic_out_ch = runDeepSomaticSingularity.out.f
        }
        if (deepsomatic_containerization == "docker") {
            deepsomatic_out_ch = runDeepSomaticDocker.out.f
        }

    emit:
        deepsomatic_out_ch
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =========================================================================================================
             Identify somatic small variants (SNVs and INDELs) in long-read DNA sequencing BAM files using DeepSomatic
             =========================================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run DeepSomatic.

        usage: nexus run --nf-workflow variant_calling_deepsomatic.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --deepsomatic_input_path            :   DeepSomatic input path.
            --deepsomatic_output_path           :   DeepSomatic output path.

        optional arguments:
            --deepsomatic_containerization      :   Containerization ('singularity' or 'docker'; default: 'singularity').
            --deepsomatic_model_type            :   DeepSomatic --model_type parameter value (default: 'PACBIO').
            --deepsomatic_bin_path              :   DeepSomatic bin path (default: 'run_deepsomatic').
            --deepsomatic_bin_version           :   DeepSomatic bin version (default: '1.9.0').
        """.stripIndent()
        exit 0
    }

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        deepsomatic_containerization        :   ${params.deepsomatic_containerization}
        deepsomatic_bin_path                :   ${params.deepsomatic_bin_path}
        deepsomatic_bin_version             :   ${params.deepsomatic_bin_version}
        deepsomatic_input_path              :   ${params.deepsomatic_input_path}
        deepsomatic_output_path             :   ${params.deepsomatic_output_path}
        deepsomatic_model_type              :   ${params.deepsomatic_model_type}
    """.stripIndent()

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

    VARIANT_CALLING_DEEPSOMATIC(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.deepsomatic_containerization,
        params.deepsomatic_bin_version,
        params.deepsomatic_bin_path,
        params.deepsomatic_input_path,
        params.deepsomatic_output_path,
        params.deepsomatic_model_type,
        params.output_dir
    )
}

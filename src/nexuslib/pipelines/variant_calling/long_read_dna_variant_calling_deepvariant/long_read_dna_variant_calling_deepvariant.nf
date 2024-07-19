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
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.containerization = 'singularity' // or docker
params.deepvariant_model_type = 'PACBIO'
params.deepvariant_bin_path = '/opt/deepvariant/bin/run_deepvariant'
params.deepvariant_bin_version = '1.6.1'
params.deepvariant_input_path = '/datastore/'
params.deepvariant_output_path = '/datastore/'
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         =================================================================================================
         Identify small variants (SNVs and INDELs) in long-read DNA sequencing BAM files using DeepVariant
         =================================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run DeepVariant.

    usage: nexus run --nf-workflow long_read_dna_small_variants.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --containerization                  :   Containerization ('singularity' or 'docker'; default: 'singularity').
        --deepvariant_model_type            :   DeepVariant --model_type parameter value (default: 'PACBIO').
        --deepvariant_bin_path              :   DeepVariant bin path (default: '/opt/deepvariant/bin/run_deepvariant').
        --deepvariant_bin_version           :   DeepVariant bin version (default: '1.6.0').
        --deepvariant_input_path            :   DeepVariant input path (default: '/datastore/').
        --deepvariant_output_path           :   DeepVariant output path (default: '/datastore/').
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
        deepvariant_input_path              :   ${params.deepvariant_input_path}
        deepvariant_output_path             :   ${params.deepvariant_output_path}
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
workflow LONG_READ_DNA_SMALL_VARIANTS {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        containerization
        deepvariant_bin_version
        deepvariant_bin_path
        deepvariant_input_path
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
                deepvariant_input_path,
                deepvariant_output_path,
                output_dir
            )
        }
        if (containerization == "docker") {
            runDeepVariantDocker(
                run_deepvariant_input_ch,
                reference_genome_fasta_file,
                deepvariant_model_type,
                deepvariant_bin_version,
                deepvariant_bin_path,
                deepvariant_input_path,
                deepvariant_output_path,
                output_dir
            )
        }
}

workflow {
    LONG_READ_DNA_SMALL_VARIANTS(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.containerization,
        params.deepvariant_bin_version,
        params.deepvariant_bin_path,
        params.deepvariant_input_path,
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

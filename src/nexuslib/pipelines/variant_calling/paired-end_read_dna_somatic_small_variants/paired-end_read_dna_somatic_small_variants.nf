#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 * Last updated date: Oct 4, 2023
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runGatk4Mutect2TumorNormal } from '../../../../modules/gatk4'
include { runGatk4LearnReadOrientationModel } from '../../../../modules/gatk4'
include { runGatk4GetPileupSummaries } from '../../../../modules/gatk4'
include { runGatk4CalculateContamination } from '../../../../modules/gatk4'
include { runGatk4FilterMutect2Calls } from '../../../../modules/gatk4'
include { runPicardMergeVCFs } from '../../../../modules/picard'
include { runStrelka2SomaticMode } from '../../../../modules/strelka2'
include { runDeepVariant as runDeepVariantTumor } from '../../../../modules/deepvariant'
include { runDeepVariant as runDeepVariantNormal } from '../../../../modules/deepvariant'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.reference_genome_fasta_file = ''
params.gatk4 = ''
params.gatk4_mutect2_params = ''
params.gatk4_getpileupsummaries_params = ''
params.picard = ''
params.strelka2 = ''
params.strelka2_params = ''
params.singularity = ''
params.deepvariant_bin_path = ''
params.deepvariant_bin_version = ''
params.deepvariant_lib_path = ''
params.deepvariant_model_type = ''
params.chromosomes = ''
// Optional arguments

// Step 3. Print inputs and help
log.info """\
         ===================================================================================
         Identify Somatic Small DNA Variants Using Paired-end Human DNA Sequencing BAM Files
         ===================================================================================
        PURPOSE:
            This workflow is intended for pairs of tumor and matched normal samples
            ()paired-end read DNA sequencing BAM files).

        WORKFLOW:
            1.  Run GATK4-Mutect2 (tumor and normal mode).
                Run GATK4 'LearnReadOrientationModel'.
                Run GATK4 'GetPileupSummaries'.
                Run GATK4 'CalculateContamination'.
                Run GATK4 'FilterMutectCalls'.
                Run Picard 'MergeVcfs',
            2.  Run Strelka2 (somatic mode).
            3.  Run DeepVariant.

         """.stripIndent()

if (params.help) {
    log.info"""\
        usage: nextflow run paired_end_human_dna_somatic_small_variants.nf [required] [optional] [--help]

        required arguments:
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file', 'tumor_sample_id', normal_sample_id'
            --output_dir                        :   Directory to which output files will be symlinked.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --gatk4                             :   GATK4 path.
            --gatk4_mutect2_params              :   GATK4-Mutect2 parameters (e.g. "--germline-resource /<path>/af-only-gnomad.hg38.vcf --panel-of-normals /<path>/1000g_pon.hg38.vcf.gz").
            --gatk4_getpileupsummaries_params   :   GATK4-GetPileupSummaries parameters (e.g. "-V /<path>/small_exac_common_3.hg38.vcf -L /<path>/small_exac_common_3.hg38.vcf").
            --picard                            :   Picard path.
            --strelka2                          :   Strelka2 path (bin path).
            --strelka2_params                   :   Strelka2 parameters.
            --singularity                       :   Singularity path.
            --deepvariant_bin_path              :   DeepVariant bin path.
            --deepvariant_bin_version           :   DeepVariant bin version.
            --deepvariant_lib_path              :   DeepVariant lib path.
            --deepvariant_model_type            :   DeepVariant --model_type parameter value.
            --chromosomes                       :   Chromosomes to parallelize using GATK4 (separated by comma; e.g. 'chr1,chr2,chr3').

        optional arguments:
    """.stripIndent()
    exit 0
} else {
    log.info"""\
            samples_tsv_file                    :   ${params.samples_tsv_file}
            output_dir                          :   ${params.output_dir}
            reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
            gatk4                               :   ${params.gatk4}
            gatk4_mutect2_params                :   ${params.gatk4_mutect2_params}
            gatk4_getpileupsummaries_params     :   ${params.gatk4_getpileupsummaries_params}
            picard                              :   ${params.picard}
            strelka2                            :   ${params.strelka2}
            strelka2_params                     :   ${params.strelka2_params}
            singularity                         :   ${params.singularity}
            deepvariant_bin_path                :   ${params.deepvariant_bin_path}
            deepvariant_bin_version             :   ${params.deepvariant_bin_version}
            deepvariant_lib_path                :   ${params.deepvariant_lib_path}
            deepvariant_model_type              :   ${params.deepvariant_model_type}
            chromosomes                         :   ${params.chromosomes}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.tumor_bam_file}",
        "${row.tumor_bam_bai_file}",
        "${row.normal_bam_file}",
        "${row.normal_bam_bai_file}",
        "${row.tumor_sample_id}",
        "${row.normal_sample_id}") }
    .set { input_bam_files_ch }

chromosomes_count = params.chromosomes.split(",").size()

// Step 5. Workflow
workflow PAIRED_END_HUMAN_DNA_SOMATIC_SMALL_VARIANTS {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(tumor_sample_id), val(normal_sample_id)]
        output_dir
        reference_genome_fasta_file
        gatk4
        gatk4_mutect2_params
        gatk4_getpileupsummaries_params
        picard
        strelka2
        strelka2_params
        singularity
        deepvariant_bin_path
        deepvariant_bin_version
        deepvariant_lib_path
        deepvariant_model_type
        chromosomes
        chromosomes_count

    main:
        chromosomes_ch = Channel
                        .value(chromosomes.tokenize(','))
                        .flatten()
        run_gatk4_mutect2_input_ch = input_bam_files_ch.combine(chromosomes_ch)
        run_strelka2_input_ch = input_bam_files_ch
        input_bam_files_ch
            .map{ [it[0], it[1], it[2]] }
            .set { run_deepvariant_tumor_input_ch }
        input_bam_files_ch
            .map{ [it[0], it[3], it[4]] }
            .set { run_deepvariant_normal_input_ch }

        // GATK4-Mutect2
        input_bam_files_ch
            .map{ [it[0], it[1], it[2]] }
            .set{ run_gatk4_getpileupsummaries_input_ch }
        runGatk4GetPileupSummaries(
            run_gatk4_getpileupsummaries_input_ch,
            gatk4,
            gatk4_getpileupsummaries_params
        )
        runGatk4CalculateContamination(
            runGatk4GetPileupSummaries.out.f,
            gatk4
        )
        runGatk4Mutect2TumorNormal(
            run_gatk4_mutect2_input_ch,
            reference_genome_fasta_file,
            gatk4,
            gatk4_mutect2_params
        )
        runGatk4Mutect2TumorNormal.out.f
            .groupTuple(by: [0], size: chromosomes_count)
            .set{ run_gatk4_mutect2_output_ch }
        run_gatk4_mutect2_output_ch
            .map{ [it[0], it[1], it[5]] }
            .transpose()
            .set{ run_gatk4_learn_read_orientation_model_input_ch }
        runGatk4LearnReadOrientationModel(
            run_gatk4_learn_read_orientation_model_input_ch,
            gatk4
        )
        runGatk4LearnReadOrientationModel.out.f
            .groupTuple(by: [0], size: chromosomes_count)
            .transpose()
            .set{ run_gatk4_learn_read_orientation_model_output_ch }
        run_gatk4_mutect2_output_ch
            .join(runGatk4CalculateContamination.out.f)
            .transpose()
            .join(run_gatk4_learn_read_orientation_model_output_ch, by: [0, 1])
            .map{ [it[0], it[1], it[2], it[3], it[4], it[6], it[7], it[8]] }
            .set { run_gatk4_filter_mutect2_calls_input_ch }
        runGatk4FilterMutect2Calls(
            run_gatk4_filter_mutect2_calls_input_ch,
            reference_genome_fasta_file,
            gatk4
        )
        runGatk4FilterMutect2Calls.out.f
            .groupTuple(by: [0], size: chromosomes_count)
            .set{ run_picard_merge_vcfs_input_ch }
        runPicardMergeVCFs(
            run_picard_merge_vcfs_input_ch,
            picard,
            "gatk4-mutect2",
            output_dir
        )

        // Strelka2
        runStrelka2SomaticMode(
            input_bam_files_ch,
            reference_genome_fasta_file,
            strelka2,
            strelka2_params,
            output_dir
        )

        // DeepVariant
        runDeepVariantTumor(
            run_deepvariant_tumor_input_ch,
            reference_genome_fasta_file,
            singularity,
            deepvariant_model_type,
            deepvariant_bin_version,
            deepvariant_bin_path,
            deepvariant_lib_path,
            output_dir
        )
        runDeepVariantNormal(
            run_deepvariant_normal_input_ch,
            reference_genome_fasta_file,
            singularity,
            deepvariant_model_type,
            deepvariant_bin_version,
            deepvariant_bin_path,
            deepvariant_lib_path,
            output_dir
        )
}

workflow {
    PAIRED_END_HUMAN_DNA_SOMATIC_SMALL_VARIANTS(
        input_bam_files_ch,
        params.output_dir,
        params.reference_genome_fasta_file,
        params.gatk4,
        params.gatk4_mutect2_params,
        params.gatk4_getpileupsummaries_params,
        params.picard,
        params.strelka2,
        params.strelka2_params,
        params.singularity,
        params.deepvariant_bin_path,
        params.deepvariant_bin_version,
        params.deepvariant_lib_path,
        params.deepvariant_model_type,
        params.chromosomes,
        chromosomes_count
    )
}

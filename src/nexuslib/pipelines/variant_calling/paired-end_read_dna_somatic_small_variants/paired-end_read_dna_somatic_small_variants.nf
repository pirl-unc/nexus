#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runGatk4Mutect2TumorNormal } from '../../modules/gatk4'
include { runGatk4Mutect2TumorNormalNonHumanSample } from '../../modules/gatk4'
include { runGatk4LearnReadOrientationModel } from '../../modules/gatk4'
include { runGatk4GetPileupSummaries } from '../../modules/gatk4'
include { runGatk4CalculateContamination } from '../../modules/gatk4'
include { runGatk4FilterMutect2Calls } from '../../modules/gatk4'
include { runGatk4FilterMutect2CallsNonHumanSample } from '../../modules/gatk4'
include { runPicardMergeVCFs } from '../../modules/picard'
include { runStrelka2SomaticMode } from '../../modules/strelka2'
include { runDeepVariantSingularity as runDeepVariantTumorSingularity } from '../../modules/deepvariant'
include { runDeepVariantSingularity as runDeepVariantNormalSingularity } from '../../modules/deepvariant'
include { runDeepVariantDocker as runDeepVariantTumorDocker } from '../../modules/deepvariant'
include { runDeepVariantDocker as runDeepVariantNormalDocker } from '../../modules/deepvariant'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.deepvariant_input_path = ''
params.deepvariant_output_path = ''
params.output_dir = ''

// Optional arguments
params.is_human = true
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.python2 = 'python2'
params.gatk4 = 'gatk4'
params.gatk4_mutect2_params = ''
params.gatk4_getpileupsummaries_params = ''
params.gatk4_chromosomes = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM'
params.picard = 'picard'
params.strelka2 = 'strelka2'
params.strelka2_params = ' '
params.containerization = 'singularity'
params.deepvariant_bin_path = '/opt/deepvariant/bin/run_deepvariant'
params.deepvariant_bin_version = '1.6.0'
params.deepvariant_model_type = 'WGS'
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         =============================================================================
         Identify somatic small DNA variants using paired-end DNA sequencing bam files
         =============================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        human:
            1.  Run gatk4-mutect2 (tumor and normal mode).
                Run gatk4 LearnReadOrientationModel.
                Run gatk4 GetPileupSummaries.
                Run gatk4 CalculateContamination.
                Run gatk4 FilterMutectCalls.
                Run picard MergeVcfs,
            2.  Run strelka2 (somatic mode).
            3.  Run deepvariant.
        non-human:
            1.  Run gatk4 mutect2 (tumor and normal mode).
                Run gatk4 FilterMutectCalls.
                Run picard MergeVcfs,
            2.  Run strelka2 (somatic mode).
            3.  Run deepvariant.

    usage: nexus run --nf-workflow paired_end_dna_somatic_small_variants.nf [required] [optional] [--help]

    required arguments:
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file', 'tumor_sample_id', normal_sample_id'
        --deepvariant_input_path            :   DeepVariant input path (e.g. /path/to/input/).
        --deepvariant_output_path           :   DeepVariant output path (e.g. /path/to/output/).
        --output_dir                        :   Directory to which output files will be symlinked.

    optional arguments:
        --is_human                          :   true if the samples are human. false otherwise (default: true).
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --python2                           :   python2 path (default: python2).
        --gatk4                             :   gatk4 path (default: gatk4).
        --gatk4_mutect2_params              :   gatk4-mutect2 parameters (e.g. '"--germline-resource /path/af-only-gnomad.hg38.vcf --panel-of-normals /path/1000g_pon.hg38.vcf.gz "').
                                                Note that the parameters need to be wrapped in quotes and a space at the end of the string is necessary.
        --gatk4_getpileupsummaries_params   :   gatk4 GetPileupSummaries parameters (e.g. '"-V /path/small_exac_common_3.hg38.vcf -L /path/small_exac_common_3.hg38.vcf "').
                                                Note that the parameters need to be wrapped in quotes and a space at the end of the string is necessary.
        --gatk4_chromosomes                 :   gatk4 chromosomes to parallelize (separated by comma; e.g. 'chr1,chr2,chr3').
        --picard                            :   Picard path (default: picard).
        --strelka2                          :   Strelka2 bin path (default: ).
        --strelka2_params                   :   Strelka2 parameters (default: ' ').
                                                Note that the parameters need to be wrapped in quotes and a space at the end of the string is necessary.
        --containerization                  :   Containerization ('singularity' or 'docker'; default: 'singularity').
        --deepvariant_bin_path              :   DeepVariant bin path (e.g. '/opt/deepvariant/bin/run_deepvariant').
        --deepvariant_bin_version           :   DeepVariant bin version (default: '1.6.0').
        --deepvariant_model_type            :   DeepVariant --model_type parameter value (default: 'WGS').
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
            samples_tsv_file                    :   ${params.samples_tsv_file}
            output_dir                          :   ${params.output_dir}
            is_human                            :   ${params.is_human}
            reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
            python2                             :   ${params.python2}
            gatk4                               :   ${params.gatk4}
            gatk4_mutect2_params                :   ${params.gatk4_mutect2_params}
            gatk4_getpileupsummaries_params     :   ${params.gatk4_getpileupsummaries_params}
            gatk4_chromosomes                   :   ${params.gatk4_chromosomes}
            picard                              :   ${params.picard}
            strelka2                            :   ${params.strelka2}
            strelka2_params                     :   ${params.strelka2_params}
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
        "${row.tumor_bam_file}",
        "${row.tumor_bam_bai_file}",
        "${row.normal_bam_file}",
        "${row.normal_bam_bai_file}",
        "${row.tumor_sample_id}",
        "${row.normal_sample_id}") }
    .set { input_bam_files_ch }

gatk4_chromosomes_count = params.gatk4_chromosomes.split(",").size()

// Step 5. Workflows
workflow PAIRED_END_READ_DNA_SOMATIC_SMALL_VARIANTS {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(tumor_sample_id), val(normal_sample_id)]
        output_dir
        is_human
        reference_genome_fasta_file
        python2
        gatk4
        gatk4_mutect2_params
        gatk4_getpileupsummaries_params
        gatk4_chromosomes
        gatk4_chromosomes_count
        picard
        strelka2
        strelka2_params
        containerization
        deepvariant_bin_path
        deepvariant_bin_version
        deepvariant_input_path
        deepvariant_output_path
        deepvariant_model_type

    main:
        gatk4_chromosomes_ch = Channel
                                .value(gatk4_chromosomes.tokenize(','))
                                .flatten()
        run_gatk4_mutect2_input_ch = input_bam_files_ch.combine(gatk4_chromosomes_ch)
        run_strelka2_input_ch = input_bam_files_ch
        input_bam_files_ch
            .map{ [it[5], it[1], it[2]] }
            .set { run_deepvariant_tumor_input_ch }
        input_bam_files_ch
            .map{ [it[6], it[3], it[4]] }
            .set { run_deepvariant_normal_input_ch }

        // GATK4-Mutect2
        if (is_human) {
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
                .groupTuple(by: [0], size: gatk4_chromosomes_count)
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
                .groupTuple(by: [0], size: gatk4_chromosomes_count)
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
                .groupTuple(by: [0], size: gatk4_chromosomes_count)
                .set{ run_picard_merge_vcfs_input_ch }
            runPicardMergeVCFs(
                run_picard_merge_vcfs_input_ch,
                picard,
                "gatk4-mutect2",
                output_dir
            )
        } else {
            // GATK4-Mutect2
            runGatk4Mutect2TumorNormalNonHumanSample(
                run_gatk4_mutect2_input_ch,
                reference_genome_fasta_file,
                gatk4,
                gatk4_mutect2_params
            )
            runGatk4Mutect2TumorNormalNonHumanSample.out.f
                .groupTuple(by: [0], size: gatk4_chromosomes_count)
                .set{ run_gatk4_mutect2_output_ch }
            run_gatk4_mutect2_output_ch
                .map{ [it[0], it[1], it[2], it[3], it[4]] }
                .set{ run_gatk4_filter_input_ch }
            runGatk4FilterMutect2CallsNonHumanSample(
                run_gatk4_filter_input_ch,
                gatk4,
                reference_genome_fasta_file
            )
            runGatk4FilterMutect2CallsNonHumanSample.out.f
                .groupTuple(by: [0], size: gatk4_chromosomes_count)
                .set{ run_picard_merge_vcfs_input_ch }
            runPicardMergeVCFs(
                run_picard_merge_vcfs_input_ch,
                picard,
                "gatk4-mutect2",
                output_dir
            )
        }

        // Strelka2
        runStrelka2SomaticMode(
            input_bam_files_ch,
            reference_genome_fasta_file,
            python2,
            strelka2,
            strelka2_params,
            output_dir
        )

        // DeepVariant
        if (containerization == "singularity") {
            runDeepVariantTumorSingularity(
                run_deepvariant_tumor_input_ch,
                reference_genome_fasta_file,
                containerization,
                deepvariant_model_type,
                deepvariant_bin_version,
                deepvariant_bin_path,
                deepvariant_input_path,
                output_dir
            )
            runDeepVariantNormalSingularity(
                run_deepvariant_normal_input_ch,
                reference_genome_fasta_file,
                containerization,
                deepvariant_model_type,
                deepvariant_bin_version,
                deepvariant_bin_path,
                deepvariant_input_path,
                output_dir
            )
        }
        if (containerization == "docker") {
            runDeepVariantTumorDocker(
                run_deepvariant_tumor_input_ch,
                reference_genome_fasta_file,
                containerization,
                deepvariant_model_type,
                deepvariant_bin_version,
                deepvariant_bin_path,
                deepvariant_input_path,
                deepvariant_output_path,
                output_dir
            )
            runDeepVariantNormalDocker(
                run_deepvariant_normal_input_ch,
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
}

workflow {
    PAIRED_END_READ_DNA_SOMATIC_SMALL_VARIANTS(
        input_bam_files_ch,
        params.output_dir,
        params.is_human,
        params.reference_genome_fasta_file,
        params.python2,
        params.gatk4,
        params.gatk4_mutect2_params,
        params.gatk4_getpileupsummaries_params,
        params.gatk4_chromosomes,
        gatk4_chromosomes_count,
        params.picard,
        params.strelka2,
        params.strelka2_params,
        params.containerization,
        params.deepvariant_bin_path,
        params.deepvariant_bin_version,
        params.deepvariant_input_path,
        params.deepvariant_output_path,
        params.deepvariant_model_type
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

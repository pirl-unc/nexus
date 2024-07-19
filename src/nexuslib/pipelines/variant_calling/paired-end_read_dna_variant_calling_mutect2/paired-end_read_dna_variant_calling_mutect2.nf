#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runGatk4Mutect2TumorNormal } from '../../modules/gatk4'
include { runGatk4LearnReadOrientationModel } from '../../modules/gatk4'
include { runGatk4GetPileupSummaries } from '../../modules/gatk4'
include { runGatk4CalculateContamination } from '../../modules/gatk4'
include { runGatk4FilterMutect2Calls } from '../../modules/gatk4'
include { runGatk4FilterMutect2CallsNonHumanSample } from '../../modules/gatk4'
include { runPicardMergeVCFs } from '../../modules/picard'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''

// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.reference_genome_fasta_dict_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.dict'
params.is_human = true
params.params_gatk4mutect2 = '--germline-resource /datastore/lbcfs/collaborations/pirl/seqdata/references/af-only-gnomad.hg38.vcf --panel-of-normals /datastore/lbcfs/collaborations/pirl/seqdata/references/1000g_pon.hg38.vcf'
params.params_gatk4getpileupsummaries = '-V /datastore/lbcfs/collaborations/pirl/seqdata/references/small_exac_common_3.hg38.vcf -L /datastore/lbcfs/collaborations/pirl/seqdata/references/small_exac_common_3.hg38.vcf'
params.chromosomes = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM'
params.delete_work_dir = false

if (params.params_gatk4mutect2 == true) {
    params_gatk4mutect2 = ''
} else {
    params_gatk4mutect2 = params.params_gatk4mutect2
}

if (params.params_gatk4getpileupsummaries == true) {
    params_gatk4getpileupsummaries = ''
} else {
    params_gatk4getpileupsummaries = params.params_gatk4getpileupsummaries
}

// Step 3. Print inputs and help
log.info """\
         =========================================================================================
         Identify somatic variants in paired-end read DNA sequencing BAM files using GATK4-Mutect2
         =========================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        human:
            1.  Run GATK4 Mutect2 (tumor and normal mode).
                Run GATK4 LearnReadOrientationModel.
                Run GATK4 GetPileupSummaries.
                Run GATK4 CalculateContamination.
                Run GATK4 FilterMutectCalls.
                Run Picard MergeVcfs,
        non-human:
            1.  Run GATK4 Mutect2 (tumor and normal mode).
                Run GATK4 FilterMutectCalls.
                Run Picard MergeVcfs,

    usage: nexus run --nf-workflow paired-end_read_dna_variant_calling_mutect2.nf [required] [optional] [--help]

    required arguments:
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file', 'tumor_sample_id', normal_sample_id'
        --output_dir                        :   Directory to which output files will be symlinked.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
        --reference_genome_fasta_dict_file  :   Reference genome DICT file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.dict).
        --is_human                          :   true if the samples are human. false otherwise (default: true).
        --params_gatk4mutect2               :   GATK4 Mutect2 parameters (default:
                                                '"--germline-resource /datastore/lbcfs/collaborations/pirl/seqdata/references/af-only-gnomad.hg38.vcf
                                                  --panel-of-normals /datastore/lbcfs/collaborations/pirl/seqdata/references/1000g_pon.hg38.vcf"').
                                                Note that the parameters need to be wrapped in quotes.
        --params_gatk4getpileupsummaries    :   GATK4 GetPileupSummaries parameters (default:
                                                '"-V /datastore/lbcfs/collaborations/pirl/seqdata/references/small_exac_common_3.hg38.vcf
                                                  -L /datastore/lbcfs/collaborations/pirl/seqdata/references/small_exac_common_3.hg38.vcf"').
                                                Note that the parameters need to be wrapped in quotes.
        --chromosomes                       :   Chromosomes to parallelize (default: 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM').
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        is_human                            :   ${params.is_human}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genome_fasta_fai_file     :   ${params.reference_genome_fasta_fai_file}
        reference_genome_fasta_dict_file    :   ${params.reference_genome_fasta_dict_file}
        params_gatk4mutect2                 :   ${params_gatk4mutect2}
        params_gatk4getpileupsummaries      :   ${params_gatk4getpileupsummaries}
        chromosomes                         :   ${params.chromosomes}
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
        "${row.normal_sample_id}") }
    .set { input_bam_files_ch }

// Step 5. Workflows
workflow PAIRED_END_READ_DNA_VARIANT_CALLING_MUTECT2 {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(tumor_sample_id), val(normal_sample_id)]
        output_dir
        reference_genome_fasta_file
        reference_genome_fasta_fai_file
        reference_genome_fasta_dict_file
        is_human
        params_gatk4mutect2
        params_gatk4getpileupsummaries
        chromosomes

    main:
        chromosomes_count = chromosomes.split(",").size()
        chromosomes_ch = Channel
            .value(chromosomes.tokenize(','))
            .flatten()
        run_gatk4_mutect2_input_ch = input_bam_files_ch.combine(chromosomes_ch)
        if (is_human) {
            input_bam_files_ch
                .map{ [it[0], it[1], it[2]] }
                .set{ run_gatk4_getpileupsummaries_input_ch }
            runGatk4GetPileupSummaries(
                run_gatk4_getpileupsummaries_input_ch,
                params_gatk4getpileupsummaries
            )
            runGatk4CalculateContamination(
                runGatk4GetPileupSummaries.out.f
            )
            runGatk4Mutect2TumorNormal(
                run_gatk4_mutect2_input_ch,
                reference_genome_fasta_file,
                reference_genome_fasta_fai_file,
                reference_genome_fasta_dict_file,
                params_gatk4mutect2
            )
            runGatk4Mutect2TumorNormal.out.f
                .groupTuple(by: [0], size: chromosomes_count)
                .set{ run_gatk4_mutect2_output_ch }
            run_gatk4_mutect2_output_ch
                .map{ [it[0], it[1], it[5]] }
                .transpose()
                .set{ run_gatk4_learn_read_orientation_model_input_ch }
            runGatk4LearnReadOrientationModel(
                run_gatk4_learn_read_orientation_model_input_ch
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
                reference_genome_fasta_fai_file,
                reference_genome_fasta_dict_file
            )
            runGatk4FilterMutect2Calls.out.f
                .groupTuple(by: [0], size: chromosomes_count)
                .set{ run_picard_merge_vcfs_input_ch }
            runPicardMergeVCFs(
                run_picard_merge_vcfs_input_ch,
                "gatk4-mutect2",
                output_dir
            )
        } else {
            runGatk4Mutect2TumorNormal(
                run_gatk4_mutect2_input_ch,
                reference_genome_fasta_file,
                reference_genome_fasta_fai_file,
                reference_genome_fasta_dict_file,
                params_gatk4mutect2
            )
            runGatk4Mutect2TumorNormal.out.f
                .groupTuple(by: [0], size: chromosomes_count)
                .set{ run_gatk4_mutect2_output_ch }
            run_gatk4_mutect2_output_ch
                .map{ [it[0], it[1], it[2], it[3], it[4]] }
                .transpose()
                .set{ run_gatk4_filter_input_ch }
            runGatk4FilterMutect2CallsNonHumanSample(
                run_gatk4_filter_input_ch,
                reference_genome_fasta_file,
                reference_genome_fasta_fai_file,
                reference_genome_fasta_dict_file
            )
            runGatk4FilterMutect2CallsNonHumanSample.out.f
                .groupTuple(by: [0], size: chromosomes_count)
                .set{ run_picard_merge_vcfs_input_ch }
            runPicardMergeVCFs(
                run_picard_merge_vcfs_input_ch,
                "gatk4-mutect2",
                output_dir
            )
        }
}

workflow {
    PAIRED_END_READ_DNA_VARIANT_CALLING_MUTECT2(
        input_bam_files_ch,
        params.output_dir,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_fai_file,
        params.reference_genome_fasta_dict_file,
        params.is_human,
        params_gatk4mutect2,
        params_gatk4getpileupsummaries,
        params.chromosomes
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

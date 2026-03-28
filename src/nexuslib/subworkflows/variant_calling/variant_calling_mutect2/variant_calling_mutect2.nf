#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }                                from '../../../tools/samtools'
include { decompressFile as decompressFasta }               from '../../../tools/utils'
include { runGatk4IndexFeature as runGatk4IndexFeature1 }   from '../../../tools/gatk4'
include { runGatk4IndexFeature as runGatk4IndexFeature2 }   from '../../../tools/gatk4'
include { runGatk4IndexFeature as runGatk4IndexFeature3 }   from '../../../tools/gatk4'
include { runGatk4CreateSequenceDictionary }                from '../../../tools/gatk4'
include { runGatk4Mutect2TumorNormal }                      from '../../../tools/gatk4'
include { runGatk4LearnReadOrientationModel }               from '../../../tools/gatk4'
include { runGatk4GetPileupSummaries }                      from '../../../tools/gatk4'
include { runGatk4CalculateContamination }                  from '../../../tools/gatk4'
include { runGatk4FilterMutect2Calls }                      from '../../../tools/gatk4'
include { runGatk4FilterMutect2CallsNonHumanSample }        from '../../../tools/gatk4'
include { runPicardMergeVCFs }                              from '../../../tools/picard'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help = ''

// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.reference_genome_fasta_file              = ''
params.mutect2_germline_resource_vcf_file       = ''
params.mutect2_panel_of_normals_vcf_file        = ''
params.getpileupsummaries_variant_vcf_file      = ''

// Optional arguments
params.params_gatk4mutect2                      = ''
params.params_gatk4getpileupsummaries           = ''
params.chromosomes                              = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM'

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

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         =========================================================================================
         Identify somatic variants in paired-end read DNA sequencing BAM files using GATK4-Mutect2
         =========================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run GATK4 Mutect2 (tumor and normal mode).
        2. Run GATK4 LearnReadOrientationModel (if --getpileupsummaries_variant_vcf_file is not empty).
        3. Run GATK4 GetPileupSummaries (if --getpileupsummaries_variant_vcf_file is not empty).
        4. Run GATK4 CalculateContamination (if --getpileupsummaries_variant_vcf_file is not empty).
        5. Run GATK4 FilterMutectCalls.
        6. Run Picard MergeVcfs,

    usage: nexus run --nf-workflow variant_calling_mutect2.nf [required] [optional] [--help]

    required arguments:
        --samples_tsv_file                      :   TSV file with the following columns:
                                                    'sample_id',
                                                    'tumor_bam_file',
                                                    'tumor_bam_bai_file',
                                                    'normal_bam_file',
                                                    'normal_bam_bai_file',
                                                    'normal_sample_id'.
        --output_dir                            :   Directory to which output files will be symlinked.
        --reference_genome_fasta_file           :   Reference genome FASTA file.
        --mutect2_germline_resource_vcf_file    :   Germline resource VCF file. This VCF file will be supplied to gatk Mutect2 --germline-resource parameter.
        --mutect2_panel_of_normals_vcf_file     :   Panel of normals VCF file. This VCF file will be supplied to gatk Mutect2 --panel-of-normals parameter.
        --getpileupsummaries_variant_vcf_file   :   GetPileupSummaries variant VCF file.

    optional arguments:
        --params_gatk4mutect2                   :   GATK4 Mutect2 parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
        --params_gatk4getpileupsummaries        :   GATK4 GetPileupSummaries parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
        --chromosomes                           :   Chromosomes to parallelize
                                                    (default: 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM').
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                        :   ${params.samples_tsv_file}
        output_dir                              :   ${params.output_dir}
        reference_genome_fasta_file             :   ${params.reference_genome_fasta_file}
        mutect2_germline_resource_vcf_file      :   ${params.mutect2_germline_resource_vcf_file}
        mutect2_panel_of_normals_vcf_file       :   ${params.mutect2_panel_of_normals_vcf_file}
        getpileupsummaries_variant_vcf_file     :   ${params.getpileupsummaries_variant_vcf_file}
        params_gatk4mutect2                     :   ${params_gatk4mutect2}
        params_gatk4getpileupsummaries          :   ${params_gatk4getpileupsummaries}
        chromosomes                             :   ${params.chromosomes}
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
        "${row.tumor_bam_file}",
        "${row.tumor_bam_bai_file}",
        "${row.normal_bam_file}",
        "${row.normal_bam_bai_file}",
        "${row.normal_sample_id}") }
    .set { input_bam_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_MUTECT2 {
    take:
        input_bam_files_ch          // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(normal_sample_id)]
        reference_genome_fasta_file
        mutect2_germline_resource_vcf_file
        mutect2_panel_of_normals_vcf_file
        getpileupsummaries_variant_vcf_file
        params_gatk4mutect2
        params_gatk4getpileupsummaries
        chromosomes
        output_dir

    main:
        // Step 1. Create input channels
        chromosomes_list            = chromosomes.tokenize(",")
        chromosomes_count           = chromosomes_list.size()
        run_gatk4_mutect2_input_ch  = input_bam_files_ch.combine(Channel.value(chromosomes_list).flatten())

        // Step 2. Decompress and index reference genome FASTA file
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidx(decompressFasta.out.f)
        runGatk4CreateSequenceDictionary(runSamtoolsFaidx.out.fasta)
        fasta_file                  = runSamtoolsFaidx.out.fasta
        fasta_fai_file              = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file              = runSamtoolsFaidx.out.gzi_file
        fasta_dict_file             = runGatk4CreateSequenceDictionary.out.f

        // Step 3. Index reference VCF files
        if (mutect2_germline_resource_vcf_file == '' || mutect2_germline_resource_vcf_file == true) {
            mutect2_germline_resource_vcf_file = file('NO_GERMLINE_FILE')
            mutect2_germline_resource_vcf_index_file = file('NO_GERMLINE_INDEX_FILE')
        } else {
            runGatk4IndexFeature1(mutect2_germline_resource_vcf_file)
            mutect2_germline_resource_vcf_file = file(mutect2_germline_resource_vcf_file)
            mutect2_germline_resource_vcf_index_file = runGatk4IndexFeature1.out.idx.mix(runGatk4IndexFeature1.out.tbi).collect()
        }
        if (mutect2_panel_of_normals_vcf_file == '' || mutect2_panel_of_normals_vcf_file == true) {
            mutect2_panel_of_normals_vcf_file = file('NO_PON_FILE')
            mutect2_panel_of_normals_vcf_index_file = file('NO_PON_INDEX_FILE')
        } else {
            runGatk4IndexFeature2(mutect2_panel_of_normals_vcf_file)
            mutect2_panel_of_normals_vcf_file = file(mutect2_panel_of_normals_vcf_file)
            mutect2_panel_of_normals_vcf_index_file = runGatk4IndexFeature2.out.idx.mix(runGatk4IndexFeature2.out.tbi).collect()
        }

        // Step 3. Run Mutect2
        if (getpileupsummaries_variant_vcf_file == '' || getpileupsummaries_variant_vcf_file == true) {
            // Run GATK4 Mutect2
            runGatk4Mutect2TumorNormal(
                run_gatk4_mutect2_input_ch,
                fasta_file,
                fasta_fai_file,
                fasta_gzi_file,
                fasta_dict_file,
                mutect2_germline_resource_vcf_file,
                mutect2_germline_resource_vcf_index_file,
                mutect2_panel_of_normals_vcf_file,
                mutect2_panel_of_normals_vcf_index_file,
                params_gatk4mutect2
            )

            // Run GATK4 FilterMutectCalls
            runGatk4Mutect2TumorNormal.out.f
                .groupTuple(by: [0], size: chromosomes_count)
                .set{ run_gatk4_mutect2_output_ch }
            run_gatk4_mutect2_output_ch
                .map{ [it[0], it[1], it[2], it[3], it[4]] }
                .transpose()
                .set{ run_gatk4_filter_input_ch }
            runGatk4FilterMutect2CallsNonHumanSample(
                run_gatk4_filter_input_ch,
                fasta_file,
                fasta_fai_file,
                fasta_gzi_file,
                fasta_dict_file
            )

            // Run Picard MergeVcfs
            runGatk4FilterMutect2CallsNonHumanSample.out.f
                .groupTuple(by: [0], size: chromosomes_count)
                .set{ run_picard_merge_vcfs_input_ch }
            runPicardMergeVCFs(
                run_picard_merge_vcfs_input_ch,
                "gatk4-mutect2",
                output_dir
            )
        } else {
            input_bam_files_ch
                .map{ [it[0], it[1], it[2]] }
                .set{ run_gatk4_getpileupsummaries_input_ch }

            // Run GATK4 GetPileupSummaries
            runGatk4IndexFeature3(getpileupsummaries_variant_vcf_file)
            getpileupsummaries_variant_vcf_file = file(getpileupsummaries_variant_vcf_file)
            getpileupsummaries_variant_vcf_index_file = runGatk4IndexFeature3.out.idx.mix(runGatk4IndexFeature3.out.tbi).collect()
            runGatk4GetPileupSummaries(
                run_gatk4_getpileupsummaries_input_ch,
                getpileupsummaries_variant_vcf_file,
                getpileupsummaries_variant_vcf_index_file,
                params_gatk4getpileupsummaries
            )

            // Run GATK4 CalculateContamination
            runGatk4CalculateContamination(runGatk4GetPileupSummaries.out.f)

            // Run GATK4 Mutect2
            runGatk4Mutect2TumorNormal(
                run_gatk4_mutect2_input_ch,
                fasta_file,
                fasta_fai_file,
                fasta_gzi_file,
                fasta_dict_file,
                mutect2_germline_resource_vcf_file,
                mutect2_germline_resource_vcf_index_file,
                mutect2_panel_of_normals_vcf_file,
                mutect2_panel_of_normals_vcf_index_file,
                params_gatk4mutect2
            )

            // Run GATK4 LearnReadOrientationModel
            runGatk4Mutect2TumorNormal.out.f
                .groupTuple(by: [0], size: chromosomes_count)
                .set{ run_gatk4_mutect2_output_ch }
            run_gatk4_mutect2_output_ch
                .map{ [it[0], it[1], it[5]] }
                .transpose()
                .set{ run_gatk4_learn_read_orientation_model_input_ch }
            runGatk4LearnReadOrientationModel(run_gatk4_learn_read_orientation_model_input_ch)

            // Run GATK4 FilterMutectCalls
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
                fasta_file,
                fasta_fai_file,
                fasta_gzi_file,
                fasta_dict_file
            )
            runGatk4FilterMutect2Calls.out.f
                .groupTuple(by: [0], size: chromosomes_count)
                .set{ run_picard_merge_vcfs_input_ch }

            // Run Picard MergeVcfs
            runPicardMergeVCFs(
                run_picard_merge_vcfs_input_ch,
                "gatk4-mutect2",
                output_dir
            )
        }

    emit:
        runPicardMergeVCFs.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_MUTECT2(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.mutect2_germline_resource_vcf_file,
        params.mutect2_panel_of_normals_vcf_file,
        params.getpileupsummaries_variant_vcf_file,
        params_gatk4mutect2,
        params_gatk4getpileupsummaries,
        params.chromosomes,
        params.output_dir
    )
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runBwaMem2 } from '../../modules/bwa'
include { runAbra2 } from '../../modules/abra2'
include { runSamtoolsSamToBam } from '../../modules/samtools'
include { runSamtoolsMarkdup } from '../../modules/samtools'
include { runSamtoolsFixmate } from '../../modules/samtools'
include { runGatk4BaseRecalibrator } from '../../modules/gatk4'
include { runGatk4GatherBQSRReports } from '../../modules/gatk4'
include { runGatk4ApplyBQSR } from '../../modules/gatk4'
include { copyBamFile } from '../../modules/utils'

// Step 2. Input parameters
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.reference_genome_fasta_file = ''
params.abra2_targets_bed_file = ''
params.bwa_mem2 = ''
params.samtools = ''
params.abra2 = ''
params.gatk4 = ''
params.gatk4_baserecalibrator_params = ''
params.chromosomes = ''
// Optional arguments
params.platform_tag = 'illumina'
params.platform_unit_tag = 'unknown'
params.library_tag = 'unknown'
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         ==========================================================
         Align Paired-end DNA Sequencing FASTQ Files Using BWA-MEM2
         ==========================================================
         PURPOSE:
            This workflow is intended for paired-end DNA sequencing FASTQ files.

         WORKFLOW:
            1. Align paired-end reads to a reference genome using bwa-mem2.
            2. Sort SAM to BAM file.
            3. Perform local realignment using abra2.
            4. Add mate score tags using samtools.
            5. Mark PCR duplicates using samtools.
            6. Calculate base recalibration scores using GATK4.
            7. Apply base recalibration scores using GATK4.
         """.stripIndent()

if (params.help) {
    log.info"""\
        usage: nextflow run paired_end_read_dna_alignment_bwa-mem2.nf [required] [optional] [--help]

        required arguments:
            --samples_tsv_file              :   TSV file with the following columns:
                                                'sample_id', 'fastq_file_1', 'fastq_file_2'.
            --output_dir                    :   Directory to which output files will be copied.
            --reference_genome_fasta_file   :   Reference genome FASTA file.
            --bwa_mem2                      :   bwa-mem2 path.
            --samtools                      :   samtools path.
            --abra2                         :   abra2 (.jar) path.
            --abra2_targets_bed_file        :   ABRA2 targets BED file.
            --gatk4                         :   GATK4 path.
            --gatk4_baserecalibrator_params :   GATK4-BaseRecalibrator parameters (e.g. '--known-sites /<path>/dbsnp_146.hg38.vcf').
            --chromosomes                   :   Chromosomes to recalibrate using GATK4 (separated by comma; e.g. 'chr1,chr2,chr3').

        optional arguments:
            --platform_tag                  :   Platform tag (default: 'illumina').
            --platform_unit_tag             :   Platform unit tag (default: 'unknown').
            --library_tag                   :   Library tag (default: 'unknown').
            --delete_work_dir               :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
            samples_tsv_file                :   ${params.samples_tsv_file}
            output_dir                      :   ${params.output_dir}
            reference_genome_fasta_file     :   ${params.reference_genome_fasta_file}
            bwa_mem2                        :   ${params.bwa_mem2}
            samtools                        :   ${params.samtools}
            abra2                           :   ${params.abra2}
            abra2_targets_bed_file          :   ${params.abra2_targets_bed_file}
            gatk4                           :   ${params.gatk4}
            gatk4_baserecalibrator_params   :   ${params.gatk4_baserecalibrator_params}
            chromosomes                     :   ${params.chromosomes}
            platform_tag                    :   ${params.platform_tag}
            platform_unit_tag               :   ${params.platform_unit_tag}
            library_tag                     :   ${params.library_tag}
            delete_work_dir                 :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.fastq_file_1}",
        "${row.fastq_file_2}") }
    .set { input_fastq_files_ch }

chromosomes_count = params.chromosomes.split(",").size()

// Step 5. Workflow
workflow PAIRED_END_DNA_READ_ALIGNMENT_BWA_MEM2 {
    take:
        input_fastq_files_ch
        reference_genome_fasta_file
        abra2_targets_bed_file
        bwa_mem2
        samtools
        abra2
        gatk4
        gatk4_baserecalibrator_params
        platform_tag
        platform_unit_tag
        library_tag
        chromosomes
        chromosomes_count
        output_dir
    main:
        chromosomes_ch = Channel
            .value(chromosomes.tokenize(','))
            .flatten()
        run_bwa_mem2_input_ch = input_fastq_files_ch

        runBwaMem2(
            run_bwa_mem2_input_ch,
            bwa_mem2,
            samtools,
            reference_genome_fasta_file,
            platform_tag,
            platform_unit_tag,
            library_tag
        )
        runSamtoolsSamToBam(
            runBwaMem2.out.f,
            samtools
        )
        runAbra2(
            runSamtoolsSamToBam.out.f,
            abra2,
            samtools,
            reference_genome_fasta_file,
            abra2_targets_bed_file
        )
        runSamtoolsFixmate(
            runAbra2.out.f,
            samtools
        )
        runSamtoolsMarkdup(
            runSamtoolsFixmate.out.f,
            samtools
        )
        runSamtoolsMarkdup.out.f
            .combine(chromosomes_ch)
            .set{ run_gatk4_base_calibrator_input_ch }
        runGatk4BaseRecalibrator(
            run_gatk4_base_calibrator_input_ch,
            reference_genome_fasta_file,
            gatk4,
            gatk4_baserecalibrator_params
        )
        runGatk4BaseRecalibrator.out.f.set { run_gatk4_base_recalibrator_output_ch }
        runGatk4BaseRecalibrator.out.f
          .groupTuple(by: [0], size: chromosomes_count)
          .map{ [it[0], it[3]] }
          .set{ run_gatk4_gather_bqsr_reports_input_ch }
        runGatk4GatherBQSRReports(
            run_gatk4_gather_bqsr_reports_input_ch,
            gatk4
        )
        run_gatk4_base_recalibrator_output_ch
           .groupTuple(by: [0])
           .map{ [it[0], it[1][0]] }
           .join(runGatk4GatherBQSRReports.out.f)
           .set{ run_gatk4_apply_bqsr_input_ch }
        runGatk4ApplyBQSR(
            run_gatk4_apply_bqsr_input_ch,
            reference_genome_fasta_file,
            gatk4,
            samtools
        )
        runGatk4ApplyBQSR.out.f.set{ copy_bam_file_input_ch }
        copyBamFile(
            copy_bam_file_input_ch,
            output_dir
        )
    emit:
        copyBamFile.out.f
}

workflow {
    PAIRED_END_DNA_READ_ALIGNMENT_BWA_MEM2(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.abra2_targets_bed_file,
        params.bwa_mem2,
        params.samtools,
        params.abra2,
        params.gatk4,
        params.gatk4_baserecalibrator_params,
        params.platform_tag,
        params.platform_unit_tag,
        params.library_tag,
        params.chromosomes,
        chromosomes_count,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
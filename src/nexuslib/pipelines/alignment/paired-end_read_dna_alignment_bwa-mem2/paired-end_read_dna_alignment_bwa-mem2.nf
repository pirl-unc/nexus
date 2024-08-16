#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runBwaMem2 } from '../../modules/bwamem2'
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
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.reference_genome_fasta_gzi_file = ''
params.reference_genome_fasta_dict_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.dict'
params.reference_genome_fasta_0123_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.0123'
params.reference_genome_fasta_amb_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.amb'
params.reference_genome_fasta_ann_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.ann'
params.reference_genome_fasta_bwt_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.bwt.2bit.64'
params.reference_genome_fasta_pac_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.pac'
params.abra2_targets_bed_file = '/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/abra2/gencode-v41-annotation-abra2-exon-targets.bed'
params.known_sites_files = '/datastore/lbcfs/collaborations/pirl/seqdata/references/dbsnp_146.hg38.vcf,/datastore/lbcfs/collaborations/pirl/seqdata/references/1000G_phase1.snps.high_confidence.hg38.vcf,/datastore/lbcfs/collaborations/pirl/seqdata/references/Mills_and_1000G_gold_standard.indels.hg38.vcf,/datastore/lbcfs/collaborations/pirl/seqdata/references/Homo_sapiens_assembly38.known_indels.vcf'
params.chromosomes = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM'
params.platform_tag = 'illumina'
params.platform_unit_tag = 'unknown'
params.library_tag = 'unknown'
params.perform_local_indel_realignment = true
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         ==========================================================
         Align paired-end DNA sequencing fastq files using bwa-mem2
         ==========================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Align paired-end reads to a reference genome using bwa-mem2.
        2. Sort sam to bam files.
        3. (Optional) Perform local realignment using abra2.
        4. Add mate score tags using samtools.
        5. Mark PCR duplicates using samtools.
        6. Calculate base recalibration scores using gatk4.
        7. Apply base recalibration scores using gatk4.

    usage: nexus run --nf-workflow paired-end_read_dna_alignment_bwa-mem2.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'fastq_file_1', 'fastq_file_2'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
        --reference_genome_fasta_gzi_file   :   Reference genome FASTA.GZ.GZI file (default: '').
        --reference_genome_fasta_dict_file  :   Reference genome DICT file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.dict).
        --reference_genome_fasta_0123_file  :   Reference genome FASTA.0123 file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.0123).
        --reference_genome_fasta_amb_file   :   Reference genome FASTA.AMB file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.amb).
        --reference_genome_fasta_ann_file   :   Reference genome FASTA.ANN file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.ann).
        --reference_genome_fasta_bwt_file   :   Reference genome FASTA.BWT file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.bwt.2bit.64).
        --reference_genome_fasta_pac_file   :   Reference genome FASTA.PAC file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.pac).
        --abra2_targets_bed_file            :   ABRA2 targets BED file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/abra2/gencode-v41-annotation-abra2-exon-targets.bed).
        --known_sites_files                 :   GATK4 BaseRecalibrator --known-sites files
                                                (default: '/datastore/lbcfs/collaborations/pirl/seqdata/references/dbsnp_146.hg38.vcf,/datastore/lbcfs/collaborations/pirl/seqdata/references/1000G_phase1.snps.high_confidence.hg38.vcf,/datastore/lbcfs/collaborations/pirl/seqdata/references/Mills_and_1000G_gold_standard.indels.hg38.vcf,/datastore/lbcfs/collaborations/pirl/seqdata/references/Homo_sapiens_assembly38.known_indels.vcf').
                                                Note that files should be separated by commas.
        --chromosomes                       :   Chromosomes to recalibrate using GATK4 (default:
                                                'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM').
        --platform_tag                      :   Platform tag (default: 'illumina').
        --platform_unit_tag                 :   Platform unit tag (default: 'unknown').
        --library_tag                       :   Library tag (default: 'unknown').
        --perform_local_indel_realignment   :   Perform local INDEL realignment (default: true).
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genome_fasta_fai_file     :   ${params.reference_genome_fasta_fai_file}
        reference_genome_fasta_gzi_file     :   ${params.reference_genome_fasta_gzi_file}
        reference_genome_fasta_dict_file    :   ${params.reference_genome_fasta_dict_file}
        reference_genome_fasta_0123_file    :   ${params.reference_genome_fasta_0123_file}
        reference_genome_fasta_amb_file     :   ${params.reference_genome_fasta_amb_file}
        reference_genome_fasta_ann_file     :   ${params.reference_genome_fasta_ann_file}
        reference_genome_fasta_bwt_file     :   ${params.reference_genome_fasta_bwt_file}
        reference_genome_fasta_pac_file     :   ${params.reference_genome_fasta_pac_file}
        abra2_targets_bed_file              :   ${params.abra2_targets_bed_file}
        known_sites_files                   :   ${params.known_sites_files}
        chromosomes                         :   ${params.chromosomes}
        platform_tag                        :   ${params.platform_tag}
        platform_unit_tag                   :   ${params.platform_unit_tag}
        library_tag                         :   ${params.library_tag}
        perform_local_indel_realignment     :   ${params.perform_local_indel_realignment}
        delete_work_dir                     :   ${params.delete_work_dir}
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

// Step 5. Workflow
workflow PAIRED_END_DNA_READ_ALIGNMENT_BWA_MEM2 {
    take:
        input_fastq_files_ch
        reference_genome_fasta_file
        reference_genome_fasta_fai_file
        reference_genome_fasta_gzi_file
        reference_genome_fasta_dict_file
        reference_genome_fasta_0123_file
        reference_genome_fasta_amb_file
        reference_genome_fasta_ann_file
        reference_genome_fasta_bwt_file
        reference_genome_fasta_pac_file
        abra2_targets_bed_file
        known_sites_files
        platform_tag
        platform_unit_tag
        library_tag
        chromosomes
        perform_local_indel_realignment
        output_dir
    main:
        if (known_sites_files == true) {
            known_sites_files_args = ''
        } else {
            known_sites_files_list = known_sites_files.split(",").collect { "--known-sites $it" }
            known_sites_files_args = known_sites_files_list.join(' ')
            println("$known_sites_files_args")
        }

        Channel
            .value(chromosomes.tokenize(","))
            .flatten()
            .set { chromosomes_ch }
        chromosomes_count = chromosomes.split(",").size()

        run_bwa_mem2_input_ch = input_fastq_files_ch

        runBwaMem2(
            run_bwa_mem2_input_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_dict_file,
            reference_genome_fasta_0123_file,
            reference_genome_fasta_amb_file,
            reference_genome_fasta_ann_file,
            reference_genome_fasta_bwt_file,
            reference_genome_fasta_fai_file,
            reference_genome_fasta_pac_file,
            platform_tag,
            platform_unit_tag,
            library_tag
        )
        runSamtoolsSamToBam(
            runBwaMem2.out.f
        )
        if (perform_local_indel_realignment == true) {
            runAbra2(
                runSamtoolsSamToBam.out.f,
                reference_genome_fasta_file,
                reference_genome_fasta_fai_file,
                abra2_targets_bed_file
            )
            runSamtoolsFixmate(
                runAbra2.out.f
            )
        } else {
            runSamtoolsFixmate(
                runSamtoolsSamToBam.out.f
            )
        }
        runSamtoolsMarkdup(
            runSamtoolsFixmate.out.f
        )
        runSamtoolsMarkdup.out.f
            .combine(chromosomes_ch)
            .set{ run_gatk4_base_calibrator_input_ch }
        runGatk4BaseRecalibrator(
            run_gatk4_base_calibrator_input_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            reference_genome_fasta_gzi_file,
            reference_genome_fasta_dict_file,
            known_sites_files_args
        )
        runGatk4BaseRecalibrator.out.f.set { run_gatk4_base_recalibrator_output_ch }
        runGatk4BaseRecalibrator.out.f
          .groupTuple(by: [0], size: chromosomes_count)
          .map{ [it[0], it[3]] }
          .set{ run_gatk4_gather_bqsr_reports_input_ch }
        runGatk4GatherBQSRReports(
            run_gatk4_gather_bqsr_reports_input_ch
        )
        run_gatk4_base_recalibrator_output_ch
           .groupTuple(by: [0])
           .map{ [it[0], it[1][0]] }
           .join(runGatk4GatherBQSRReports.out.f)
           .set{ run_gatk4_apply_bqsr_input_ch }
        runGatk4ApplyBQSR(
            run_gatk4_apply_bqsr_input_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            reference_genome_fasta_gzi_file,
            reference_genome_fasta_dict_file
        )
        runGatk4ApplyBQSR.out.f.set{ copy_bam_file_input_ch }
        copyBamFile(
            copy_bam_file_input_ch,
            output_dir
        )
    emit:
        runGatk4ApplyBQSR.out.f
}

workflow {
    PAIRED_END_DNA_READ_ALIGNMENT_BWA_MEM2(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_fai_file,
        params.reference_genome_fasta_gzi_file,
        params.reference_genome_fasta_dict_file,
        params.reference_genome_fasta_0123_file,
        params.reference_genome_fasta_amb_file,
        params.reference_genome_fasta_ann_file,
        params.reference_genome_fasta_bwt_file,
        params.reference_genome_fasta_pac_file,
        params.abra2_targets_bed_file,
        params.known_sites_files,
        params.platform_tag,
        params.platform_unit_tag,
        params.library_tag,
        params.chromosomes,
        params.perform_local_indel_realignment,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
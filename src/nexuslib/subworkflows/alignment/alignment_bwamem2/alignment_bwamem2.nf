#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }                    from '../../../tools/samtools'
include { runBwaMem2Index }                     from '../../../tools/bwamem2'
include { runGatk4CreateSequenceDictionary }    from '../../../tools/gatk4'
include { runGatk4IndexFeature }                from '../../../tools/gatk4'
include { runBwaMem2 }                          from '../../../tools/bwamem2'
include { runAbra2 }                            from '../../../tools/abra2'
include { runSamtoolsMarkdup }                  from '../../../tools/samtools'
include { runSamtoolsFixmate }                  from '../../../tools/samtools'
include { runGatk4BaseRecalibrator }            from '../../../tools/gatk4'
include { runGatk4GatherBQSRReports }           from '../../../tools/gatk4'
include { runGatk4ApplyBQSR }                   from '../../../tools/gatk4'
include { copyBamFile }                         from '../../../tools/utils'
include { decompressFile as decompressFasta }   from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.abra2_targets_bed_file           = ''
params.known_sites_vcf_files            = ''

// Optional arguments
params.chromosomes                      = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM'
params.platform_tag                     = 'illumina'
params.platform_unit_tag                = 'unknown'
params.library_tag                      = 'unknown'
params.perform_local_indel_realignment  = true

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow INDEX_REFERENCE_FASTA {
    take:
        reference_fasta_file

    main:
        runSamtoolsFaidx(reference_fasta_file)
        runBwaMem2Index(runSamtoolsFaidx.out.fasta)
        runGatk4CreateSequenceDictionary(runSamtoolsFaidx.out.fasta)

    emit:
        fasta_file          = runSamtoolsFaidx.out.fasta
        fasta_fai_file      = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file      = runSamtoolsFaidx.out.gzi_file
        fasta_dict_file     = runGatk4CreateSequenceDictionary.out.f
        fasta_0123_file     = runBwaMem2Index.out.f.map { dir -> file("${dir}/*.0123") }
        fasta_amb_file      = runBwaMem2Index.out.f.map { dir -> file("${dir}/*.amb") }
        fasta_ann_file      = runBwaMem2Index.out.f.map { dir -> file("${dir}/*.ann") }
        fasta_bwt_file      = runBwaMem2Index.out.f.map { dir -> file("${dir}/*.bwt.2bit.64") }
        fasta_pac_file      = runBwaMem2Index.out.f.map { dir -> file("${dir}/*.pac") }
}

workflow ALIGNMENT_BWAMEM2 {
    take:
        input_fastq_files_ch
        reference_genome_fasta_file
        abra2_targets_bed_file
        known_sites_vcf_files
        platform_tag
        platform_unit_tag
        library_tag
        chromosomes
        perform_local_indel_realignment
        output_dir

    main:
        // Step 1. Create input channels
        chromosomes_list    = chromosomes.tokenize(",")
        chromosomes_count   = chromosomes_list.size()

        // Step 2. Index reference genome FASTA file
        INDEX_REFERENCE_FASTA(reference_genome_fasta_file)
        fasta_file          = INDEX_REFERENCE_FASTA.out.fasta_file
        fasta_fai_file      = INDEX_REFERENCE_FASTA.out.fasta_fai_file
        fasta_gzi_file      = INDEX_REFERENCE_FASTA.out.fasta_gzi_file
        fasta_dict_file     = INDEX_REFERENCE_FASTA.out.fasta_dict_file
        fasta_0123_file     = INDEX_REFERENCE_FASTA.out.fasta_0123_file
        fasta_amb_file      = INDEX_REFERENCE_FASTA.out.fasta_amb_file
        fasta_ann_file      = INDEX_REFERENCE_FASTA.out.fasta_ann_file
        fasta_bwt_file      = INDEX_REFERENCE_FASTA.out.fasta_bwt_file
        fasta_pac_file      = INDEX_REFERENCE_FASTA.out.fasta_pac_file

        // Step 3. Index known sites VCF files
        known_sites_vcf_paths = known_sites_vcf_files.split(",").collect { it.trim() }
        Channel.from(known_sites_vcf_paths).map { file(it) }.set { known_sites_vcfs_ch }
        runGatk4IndexFeature(known_sites_vcfs_ch)
        known_sites_files       = known_sites_vcf_paths.collect { file(it) }
        known_sites_index_files = runGatk4IndexFeature.out.idx.mix(runGatk4IndexFeature.out.tbi).collect()

        // Step 4. Run Bwa-mem2 (pipes bwa-mem2 mem | samtools view | samtools sort, then samtools index)
        runBwaMem2(
            input_fastq_files_ch,
            fasta_file,
            fasta_fai_file,
            fasta_dict_file,
            fasta_0123_file,
            fasta_amb_file,
            fasta_ann_file,
            fasta_bwt_file,
            fasta_pac_file,
            platform_tag,
            platform_unit_tag,
            library_tag
        )

        // Step 5. Perform local INDEL realignment
        if (perform_local_indel_realignment == true) {
            decompressFasta(reference_genome_fasta_file)
            runAbra2(
                runBwaMem2.out.f,
                decompressFasta.out.f,
                abra2_targets_bed_file
            )
            runSamtoolsFixmate(runAbra2.out.f)
        } else {
            runSamtoolsFixmate(runBwaMem2.out.f)
        }

        // Step 6. Perform mark duplicate
        runSamtoolsMarkdup(runSamtoolsFixmate.out.f)

        // Step 7. Recalibrate base
        runSamtoolsMarkdup.out.f
            .combine(Channel.from(chromosomes_list))
            .set{ run_gatk4_base_calibrator_input_ch }
        runGatk4BaseRecalibrator(
            run_gatk4_base_calibrator_input_ch,
            fasta_file,
            fasta_fai_file,
            fasta_gzi_file,
            fasta_dict_file,
            known_sites_files,
            known_sites_index_files
        )

        // Step 8. Gather BQSR reports
        runGatk4BaseRecalibrator.out.f.set { run_gatk4_base_recalibrator_output_ch }
        runGatk4BaseRecalibrator.out.f
          .groupTuple(by: [0], size: chromosomes_count)
          .map{ [it[0], it[3]] }
          .set{ run_gatk4_gather_bqsr_reports_input_ch }
        runGatk4GatherBQSRReports(run_gatk4_gather_bqsr_reports_input_ch)

        // Step 9. Apply BQSR
        run_gatk4_base_recalibrator_output_ch
           .groupTuple(by: [0])
           .map{ [it[0], it[1][0]] }
           .join(runGatk4GatherBQSRReports.out.f)
           .set{ run_gatk4_apply_bqsr_input_ch }
        runGatk4ApplyBQSR(
            run_gatk4_apply_bqsr_input_ch,
            fasta_file,
            fasta_fai_file,
            fasta_gzi_file,
            fasta_dict_file
        )

        // Step 10. Copy BAM files
        copyBamFile(
            runGatk4ApplyBQSR.out.f,
            output_dir
        )

    emit:
        runGatk4ApplyBQSR.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ==========================================================
             Align paired-end DNA sequencing fastq files using bwa-mem2
             ==========================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Index fasta and vcf files.
            2. Align paired-end reads to a reference genome using bwa-mem2 piped to a sorted BAM
               (bwa-mem2 mem | samtools view | samtools sort), then samtools index.
            3. Perform local realignment using abra2 (optional).
            4. Add mate score tags using samtools.
            5. Mark PCR duplicates using samtools.
            6. Calculate base recalibration scores using gatk4.
            7. Apply base recalibration scores using gatk4.

        usage: nexus run --nf-workflow alignment_bwamem2.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'fastq_file_1', 'fastq_file_2'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --abra2_targets_bed_file            :   ABRA2 targets BED file.
            --known_sites_vcf_files             :   GATK4 BaseRecalibrator --known-sites files. At least one VCF file must be supplied. Note that VCF files should be separated by commas.

        optional arguments:
            --chromosomes                       :   Chromosomes to recalibrate using GATK4 (default:
                                                    'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM').
            --platform_tag                      :   Platform tag (default: 'illumina').
            --platform_unit_tag                 :   Platform unit tag (default: 'unknown').
            --library_tag                       :   Library tag (default: 'unknown').
            --perform_local_indel_realignment   :   Perform local INDEL realignment (default: true).
        """.stripIndent()
        exit 0
    }

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        abra2_targets_bed_file              :   ${params.abra2_targets_bed_file}
        known_sites_vcf_files               :   ${params.known_sites_vcf_files}
        chromosomes                         :   ${params.chromosomes}
        platform_tag                        :   ${params.platform_tag}
        platform_unit_tag                   :   ${params.platform_unit_tag}
        library_tag                         :   ${params.library_tag}
        perform_local_indel_realignment     :   ${params.perform_local_indel_realignment}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file_1}",
            "${row.fastq_file_2}") }
        .set { input_fastq_files_ch }

    ALIGNMENT_BWAMEM2(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.abra2_targets_bed_file,
        params.known_sites_vcf_files,
        params.platform_tag,
        params.platform_unit_tag,
        params.library_tag,
        params.chromosomes,
        params.perform_local_indel_realignment,
        params.output_dir
    )
}

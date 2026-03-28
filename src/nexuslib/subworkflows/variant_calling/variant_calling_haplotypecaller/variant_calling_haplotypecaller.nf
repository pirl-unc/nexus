#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }                                from '../../../tools/samtools'
include { runGatk4CreateSequenceDictionary }                from '../../../tools/gatk4'
include { runGatk4HaplotypeCaller }                         from '../../../tools/gatk4'
include { runPicardMergeVCFs }                              from '../../../tools/picard'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help = ''

// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.reference_genome_fasta_file              = ''

// Optional arguments
params.params_gatk4haplotypecaller              = ''
params.chromosomes                              = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM'

if (params.params_gatk4haplotypecaller == true) {
    params_gatk4haplotypecaller = ''
} else {
    params_gatk4haplotypecaller = params.params_gatk4haplotypecaller
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ==================================================================================================
         Identify germline variants in paired-end read DNA sequencing BAM files using GATK4-HaplotypeCaller
         ==================================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run GATK4 HaplotypeCaller.
        2. Run Picard MergeVcfs,

    usage: nexus run --nf-workflow variant_calling_haplotypecaller.nf [required] [optional] [--help]

    required arguments:
        --samples_tsv_file                      :   TSV file with the following columns:
                                                    'sample_id',
                                                    'bam_file',
                                                    'bam_bai_file'.
        --output_dir                            :   Directory to which output files will be symlinked.
        --reference_genome_fasta_file           :   Reference genome FASTA file.

    optional arguments:
        --params_gatk4haplotypecaller           :   GATK4 HaplotypeCaller parameters (default: '""').
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
        params_gatk4haplotypecaller             :   ${params_gatk4haplotypecaller}
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
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_HAPLOTYPECALLER {
    take:
        input_bam_files_ch          // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        output_dir
        reference_genome_fasta_file
        params_gatk4haplotypecaller
        chromosomes

    main:
        // Step 1. Create input channels
        chromosomes_list            = chromosomes.tokenize(",")
        chromosomes_count           = chromosomes_list.size()
        run_gatk4_haplotypecaller_input_ch = input_bam_files_ch.combine(Channel.value(chromosomes_list).flatten())

        // Step 2. Index reference genome FASTA file
        runSamtoolsFaidx(reference_genome_fasta_file)
        runGatk4CreateSequenceDictionary(runSamtoolsFaidx.out.fasta)
        fasta_file                  = runSamtoolsFaidx.out.fasta
        fasta_fai_file              = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file              = runSamtoolsFaidx.out.gzi_file
        fasta_dict_file             = runGatk4CreateSequenceDictionary.out.f

        // Step 3. Run HaplotypeCaller
        runGatk4HaplotypeCaller(
            run_gatk4_haplotypecaller_input_ch,
            fasta_file,
            fasta_fai_file,
            fasta_gzi_file,
            fasta_dict_file,
            params_gatk4haplotypecaller,
            output_dir
        )

        // Step 4. Run Picard MergeVcfs
        runGatk4HaplotypeCaller.out.f
            .groupTuple(by: [0], size: chromosomes_count)
            .set{ run_picard_merge_vcfs_input_ch }
        runPicardMergeVCFs(
            run_picard_merge_vcfs_input_ch,
            "gatk4-haplotypecaller",
            output_dir
        )

    emit:
        runPicardMergeVCFs.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_HAPLOTYPECALLER(
        input_bam_files_ch,
        params.output_dir,
        params.reference_genome_fasta_file,
        params_gatk4haplotypecaller,
        params.chromosomes
    )
}

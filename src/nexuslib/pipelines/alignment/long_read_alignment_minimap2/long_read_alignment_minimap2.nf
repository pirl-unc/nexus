#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runMinimap2 } from '../../modules/minimap2'
include { runSamtoolsSamToBam } from '../../modules/samtools'
include { runSamtoolsCalmd } from '../../modules/samtools'
include { runSamtoolsSort } from '../../modules/samtools'
include { copyBamFile } from '../../modules/utils'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.reference_genome_fasta_file = ''
params.minimap2 = ''
params.minimap2_params = ''
params.samtools = ''
// Optional arguments
params.platform_tag = 'unknown'
params.platform_unit_tag = 'unknown'
params.library_tag = 'unknown'
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         ==================================================================
         Align Long-read (DNA or RNA) Sequencing FASTQ Files Using Minimap2
         ==================================================================
         Workflow:
            1. Align reads to a reference genome using minimap2.
            2. Convert SAM to BAM.
            3. Generate MD tags.
            4. Sort MD-tagged BAM file.
         """.stripIndent()

if (params.help) {
    log.info"""\
        workflow:
            1. Align reads (fastq.gz files) to a reference genome using minimap2.
            2. Convert sam files to bam files.
            3. Generate MD tags.
            4. Sort MD-tagged bam files.

        usage: nexus run --nf-workflow long_read_alignment_minimap2.nf [required] [optional] [--help]

        required arguments:
            -c                              :   Nextflow .config file.
            -w                              :   Nextflow work directory path.
            --samples_tsv_file              :   TSV file with the following columns:
                                                'sample_id', 'fastq_file'.
            --output_dir                    :   Directory to which output files will be copied.
            --reference_genome_fasta_file   :   Reference genome FASTA file.
            --minimap2                      :   minimap2 path.
            --minimap2_params               :   minimap2 parameters (e.g. '"--ax map-hifi --cs --eqx --Y --L "').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
            --samtools                      :   samtools path.

        optional arguments:
            --platform_tag                  :   Platform tag (default: 'unknown').
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
            minimap2                        :   ${params.minimap2}
            minimap2_params                 :   ${params.minimap2_params}
            samtools                        :   ${params.samtools}
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
        "${row.fastq_file}") }
    .set { input_fastq_files_ch }

// Step 5. Workflow
workflow LONG_READ_ALIGNMENT_MINIMAP2 {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        reference_genome_fasta_file
        minimap2
        minimap2_params
        samtools
        platform_tag
        platform_unit_tag
        library_tag
        output_dir
    main:
        run_minimap2_input_ch = input_fastq_files_ch
        runMinimap2(
            run_minimap2_input_ch,
            reference_genome_fasta_file,
            minimap2,
            minimap2_params,
            samtools,
            platform_tag,
            platform_unit_tag,
            library_tag
        )
        runSamtoolsSamToBam(
            runMinimap2.out.f,
            samtools
        )
        runSamtoolsCalmd(
            runSamtoolsSamToBam.out.f,
            samtools,
            reference_genome_fasta_file
        )
        runSamtoolsSort(
            runSamtoolsCalmd.out.f,
            samtools
        )
        runSamtoolsSort.out.f.set{ copy_bam_file_input_ch }
        copyBamFile(
            copy_bam_file_input_ch,
            output_dir
        )
    emit:
        copyBamFile.out.f
}

workflow {
    LONG_READ_ALIGNMENT_MINIMAP2(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.minimap2,
        params.minimap2_params,
        params.samtools,
        params.platform_tag,
        params.platform_unit_tag,
        params.library_tag,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
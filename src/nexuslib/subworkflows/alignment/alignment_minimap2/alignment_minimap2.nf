#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }        from '../../../tools/samtools'
include { runMinimap2 }             from '../../../tools/minimap2'
include { runSamtoolsSamToBam }     from '../../../tools/samtools'
include { runSamtoolsCalmd }        from '../../../tools/samtools'
include { runSamtoolsSort }         from '../../../tools/samtools'
include { copyBamFile }             from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir = ''
params.reference_genome_fasta_file  = ''

// Optional arguments
params.params_minimap2              = '-ax map-hifi --cs --eqx -Y -L --secondary=no'
params.platform_tag                 = 'unknown'
params.platform_unit_tag            = 'unknown'
params.library_tag                  = 'unknown'

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         =======================================================
         Align long-read (DNA or RNA) fastq files using minimap2
         =======================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Align reads (fastq.gz files) to a reference genome using minimap2.
        2. Convert sam files to bam files.
        3. Generate MD tags.
        4. Sort MD-tagged bam files.

    usage: nexus run --nf-workflow alignment_minimap2.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'fastq_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.

    optional arguments:
        --params_minimap2                   :   Minimap2 parameters (default: "-ax map-hifi --cs --eqx -Y -L --secondary=no").
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
        --platform_tag                      :   Platform tag (default: 'unknown').
        --platform_unit_tag                 :   Platform unit tag (default: 'unknown').
        --library_tag                       :   Library tag (default: 'unknown').
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_minimap2                     :   ${params.params_minimap2}
        platform_tag                        :   ${params.platform_tag}
        platform_unit_tag                   :   ${params.platform_unit_tag}
        library_tag                         :   ${params.library_tag}
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
        "${row.fastq_file}") }
    .set { input_fastq_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow ALIGNMENT_MINIMAP2 {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        reference_genome_fasta_file
        params_minimap2
        platform_tag
        platform_unit_tag
        library_tag
        output_dir

    main:
        // Step 1. Index reference genome FASTA file
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file          = runSamtoolsFaidx.out.fasta
        fasta_fai_file      = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file      = runSamtoolsFaidx.out.gzi_file

        // Step 2. Run Minimap2
        run_minimap2_input_ch = input_fastq_files_ch
        runMinimap2(
            run_minimap2_input_ch,
            fasta_file,
            fasta_fai_file,
            params_minimap2,
            platform_tag,
            platform_unit_tag,
            library_tag
        )

        // Step 3. Convert SAM to BAM
        runSamtoolsSamToBam(runMinimap2.out.f)

        // Step 4. Run Samtools calmd
        runSamtoolsCalmd(
            runSamtoolsSamToBam.out.f,
            fasta_file,
            fasta_fai_file
        )

        // Step 5. Sort the BAM files
        runSamtoolsSort(runSamtoolsCalmd.out.f)

        // Step 6. Copy BAM files
        runSamtoolsSort.out.f.set{ copy_bam_file_input_ch }
        copyBamFile(
            copy_bam_file_input_ch,
            output_dir
        )

    emit:
        runSamtoolsSort.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ALIGNMENT_MINIMAP2(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.params_minimap2,
        params.platform_tag,
        params.platform_unit_tag,
        params.library_tag,
        params.output_dir
    )
}

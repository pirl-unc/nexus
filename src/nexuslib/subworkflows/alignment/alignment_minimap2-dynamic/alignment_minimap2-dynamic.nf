#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }                    from '../../../tools/samtools'
include { runMinimap2CustomReference }          from '../../../tools/minimap2'
include { runSamtoolsSamToBamCustomReference }  from '../../../tools/samtools'
include { runSamtoolsCalmdCustomReference }     from '../../../tools/samtools'
include { runSamtoolsSortCustomReference }      from '../../../tools/samtools'
include { copyBamFile }                         from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                     = ''

// Required arguments
params.samples_tsv_file         = ''
params.output_dir               = ''

// Optional arguments
params.params_minimap2          = '-ax map-hifi --cs --eqx -Y -L --secondary=no'
params.platform_tag             = 'unknown'
params.platform_unit_tag        = 'unknown'
params.library_tag              = 'unknown'

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         =================================================================================
         Align long-read (DNA or RNA) fastq files against custom references using minimap2
         =================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Align reads (fastq.gz files) to a custom reference genome using minimap2.
        2. Convert sam files to bam files.
        3. Generate MD tags.
        4. Sort MD-tagged bam files.

    usage: nexus run --nf-workflow alignment_minimap2-dynamic.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'fastq_file', 'reference_genome_fasta_file', 'reference_genome_fasta_fai_file'.
        --output_dir                        :   Directory to which output files will be copied.

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
        "${row.fastq_file}",
        "${row.reference_genome_fasta_file}",
        "${row.reference_genome_fasta_fai_file}") }
    .set { input_fastq_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow ALIGNMENT_MINIMAP2_DYNAMIC {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file), path(reference_genome_fasta_file), path(reference_genome_fasta_fai_file)]
        params_minimap2
        platform_tag
        platform_unit_tag
        library_tag
        output_dir
    main:
        run_minimap2_input_ch = input_fastq_files_ch
        runMinimap2CustomReference(
            run_minimap2_input_ch,
            params_minimap2,
            platform_tag,
            platform_unit_tag,
            library_tag
        )
        runSamtoolsSamToBamCustomReference(runMinimap2CustomReference.out.f)
        runSamtoolsCalmdCustomReference(runSamtoolsSamToBamCustomReference.out.f)
        runSamtoolsSortCustomReference(runSamtoolsCalmdCustomReference.out.f)
        runSamtoolsSortCustomReference.out.f.set{ copy_bam_file_input_ch }
        copyBamFile(
            copy_bam_file_input_ch,
            output_dir
        )
    emit:
        copyBamFile.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ALIGNMENT_MINIMAP2_DYNAMIC(
        input_fastq_files_ch,
        params.params_minimap2,
        params.platform_tag,
        params.platform_unit_tag,
        params.library_tag,
        params.output_dir
    )
}

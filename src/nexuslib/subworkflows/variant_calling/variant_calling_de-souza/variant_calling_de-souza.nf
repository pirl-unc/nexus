#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }                        from '../../../tools/samtools'
include { runGatk4CreateSequenceDictionary }        from '../../../tools/gatk4'
include { runMinimap2 }                             from '../../../tools/minimap2'
include { runSamtoolsSamToBam }                     from '../../../tools/samtools'
include { runSamtoolsFilter }                       from '../../../tools/samtools'
include { runSamtoolsSort }                         from '../../../tools/samtools'
include { runSamtoolsIndex as runSamtoolsIndex1 }   from '../../../tools/samtools'
include { runSamtoolsIndex as runSamtoolsIndex2 }   from '../../../tools/samtools'
include { runSamtoolsIndex as runSamtoolsIndex3 }   from '../../../tools/samtools'
include { runGatk4SplitNCigarReads }                from '../../../tools/gatk4'
include { runFlagCorrection }                       from '../../../tools/de_souza'
include { runDeepVariantSingularity }               from '../../../tools/deepvariant'
include { runDeepVariantDocker }                    from '../../../tools/deepvariant'
include { copyBamFile }                             from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.deepvariant_input_path           = ''
params.deepvariant_output_path          = ''

// Optional arguments
params.params_minimap2                  = '-ax splice -uf -C5 --secondary=no'
params.params_samtools_view             = '-F 2308'
params.platform_tag                     = 'unknown'
params.platform_unit_tag                = 'unknown'
params.library_tag                      = 'unknown'
params.deepvariant_containerization     = 'singularity' // or docker
params.deepvariant_model_type           = 'PACBIO'
params.deepvariant_bin_path             = '/opt/deepvariant/bin/run_deepvariant'
params.deepvariant_bin_version          = '1.6.1'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_DE_SOUZA {
    take:
        input_fastq_files_ch             // channel: [val(sample_id), path(fastq_file)]
        reference_genome_fasta_file
        params_minimap2
        platform_tag
        platform_unit_tag
        library_tag
        params_samtools_view
        deepvariant_containerization
        deepvariant_bin_version
        deepvariant_bin_path
        deepvariant_input_path
        deepvariant_output_path
        deepvariant_model_type
        output_dir

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file      = runSamtoolsFaidx.out.fasta
        fasta_fai_file  = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file  = runSamtoolsFaidx.out.gzi_file

        runGatk4CreateSequenceDictionary(reference_genome_fasta_file)
        dict_file       = runGatk4CreateSequenceDictionary.out.f

        runMinimap2(
            input_fastq_files_ch,
            fasta_file,
            fasta_fai_file,
            params_minimap2,
            platform_tag,
            platform_unit_tag,
            library_tag
        )

        runSamtoolsSamToBam(runMinimap2.out.f)

        runSamtoolsFilter(
            runSamtoolsSamToBam.out.f,
            params_samtools_view
        )

        runSamtoolsIndex1(runSamtoolsFilter.out.f)

        runSamtoolsSort(runSamtoolsIndex1.out.f)

        run_samtools_sort_output_ch = runSamtoolsSort.out.f
        runGatk4SplitNCigarReads(
            run_samtools_sort_output_ch,
            fasta_file,
            fasta_fai_file,
            fasta_gzi_file,
            dict_file
        )

        run_gatk4_splitncigarreads_output_ch = runGatk4SplitNCigarReads.out.f
        runSamtoolsIndex2(run_gatk4_splitncigarreads_output_ch)
        run_samtools_index_output_ch = runSamtoolsIndex2.out.f
        run_flag_correction_input_ch = run_samtools_sort_output_ch.join(run_samtools_index_output_ch)
            .map { sample_id, bam_file, bam_bai_file, splitncigarreads_bam_file, splitncigarreads_bam_bai_file ->
                [sample_id, bam_file, bam_bai_file, splitncigarreads_bam_file, splitncigarreads_bam_bai_file]
            }

        runFlagCorrection(run_flag_correction_input_ch)

        runSamtoolsIndex3(runFlagCorrection.out.f)

        copy_bam_file_input_ch = runSamtoolsIndex3.out.f

        run_deepvariant_input_ch = copy_bam_file_input_ch
        copyBamFile(
            copy_bam_file_input_ch,
            output_dir
        )

        // DeepVariant
        if (deepvariant_containerization == "singularity") {
            runDeepVariantSingularity(
                run_deepvariant_input_ch,
                fasta_file,
                fasta_fai_file,
                deepvariant_containerization,
                deepvariant_model_type,
                deepvariant_bin_version,
                deepvariant_bin_path,
                deepvariant_input_path,
                deepvariant_output_path,
                output_dir
            )
        }
        if (deepvariant_containerization == "docker") {
            runDeepVariantDocker(
                run_deepvariant_input_ch,
                fasta_file,
                fasta_fai_file,
                deepvariant_model_type,
                deepvariant_bin_version,
                deepvariant_bin_path,
                deepvariant_input_path,
                deepvariant_output_path,
                output_dir
            )
        }

        de_souza_out_ch = Channel.empty()
        if (deepvariant_containerization == "singularity") {
            de_souza_out_ch = runDeepVariantSingularity.out.f
        }
        if (deepvariant_containerization == "docker") {
            de_souza_out_ch = runDeepVariantDocker.out.f
        }

    emit:
        de_souza_out_ch
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =================================================================================================================================
             Identify RNA variants in long-read RNA sequencing FASTQ files using lrRNAseqVariantCalling (de Souza et al., Genome Biology 2023)
             =================================================================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Align reads to the reference using minimap2.
            2. Filter, sort, and index using samtools.
            3. Run GATK4 SplitNCigarReads.
            4. Perform flag correction.
            5. Index using samtools.
            6. Run DeepVariant.

        usage: nexus run --nf-workflow variant_calling_de-souza.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'fastq_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --deepvariant_input_path            :   DeepVariant input path.
            --deepvariant_output_path           :   DeepVariant output path.

        optional arguments:
            --params_minimap2                   :   Minimap2 parameters (default: '"-ax splice -uf -C5 --secondary=no"').
                                                    Note that the parameters need to be wrapped in quotes.
            --params_samtools_view              :   Samtools view parameters (default: '"-F 2308"').
                                                    Note that the parameters need to be wrapped in quotes.
            --platform_tag                      :   Platform tag (default: 'unknown').
            --platform_unit_tag                 :   Platform unit tag (default: 'unknown').
            --library_tag                       :   Library tag (default: 'unknown').
            --deepvariant_containerization      :   DeepVariant containerization ('singularity' or 'docker'; default: 'singularity').
            --deepvariant_model_type            :   DeepVariant --model_type parameter value (default: 'PACBIO').
            --deepvariant_bin_path              :   DeepVariant bin path (default: '/opt/deepvariant/bin/run_deepvariant').
            --deepvariant_bin_version           :   DeepVariant bin version (default: '1.6.0').
        """.stripIndent()
        exit 0
    }

    def params_minimap2      = (params.params_minimap2 == true) ? '' : params.params_minimap2
    def params_samtools_view = (params.params_samtools_view == true) ? '' : params.params_samtools_view

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_minimap2                     :   ${params_minimap2}
        params_samtools_view                :   ${params_samtools_view}
        platform_tag                        :   ${params.platform_tag}
        platform_unit_tag                   :   ${params.platform_unit_tag}
        library_tag                         :   ${params.library_tag}
        deepvariant_containerization        :   ${params.deepvariant_containerization}
        deepvariant_model_type              :   ${params.deepvariant_model_type}
        deepvariant_bin_path                :   ${params.deepvariant_bin_path}
        deepvariant_bin_version             :   ${params.deepvariant_bin_version}
        deepvariant_input_path              :   ${params.deepvariant_input_path}
        deepvariant_output_path             :   ${params.deepvariant_output_path}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_fastq_files_ch }

    VARIANT_CALLING_DE_SOUZA(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params_minimap2,
        params.platform_tag,
        params.platform_unit_tag,
        params.library_tag,
        params_samtools_view,
        params.deepvariant_containerization,
        params.deepvariant_bin_version,
        params.deepvariant_bin_path,
        params.deepvariant_input_path,
        params.deepvariant_output_path,
        params.deepvariant_model_type,
        params.output_dir
    )
}

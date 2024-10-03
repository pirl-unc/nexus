#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runMinimap2 } from '../../modules/minimap2'
include { runSamtoolsSamToBam } from '../../modules/samtools'
include { runSamtoolsFilter } from '../../modules/samtools'
include { runSamtoolsSort } from '../../modules/samtools'
include { runSamtoolsIndex as runSamtoolsIndex1 } from '../../modules/samtools'
include { runSamtoolsIndex as runSamtoolsIndex2 } from '../../modules/samtools'
include { runSamtoolsIndex as runSamtoolsIndex3 } from '../../modules/samtools'
include { runGatk4SplitNCigarReads } from '../../modules/gatk4'
include { runFlagCorrection } from '../../modules/de_souza'
include { runDeepVariantSingularity } from '../../modules/deepvariant'
include { runDeepVariantDocker } from '../../modules/deepvariant'
include { copyBamFile } from '../../modules/utils'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.reference_genome_fasta_dict_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.dict'
params.params_minimap2 = '-ax splice -uf -C5 --secondary=no'
params.params_samtools_view = '-F 2308'
params.platform_tag = 'unknown'
params.platform_unit_tag = 'unknown'
params.library_tag = 'unknown'
params.deepvariant_containerization = 'singularity' // or docker
params.deepvariant_model_type = 'PACBIO'
params.deepvariant_bin_path = '/opt/deepvariant/bin/run_deepvariant'
params.deepvariant_bin_version = '1.6.1'
params.deepvariant_input_path = '/datastore/'
params.deepvariant_output_path = '/datastore/'
params.delete_work_dir = false

if (params.params_minimap2 == true) {
    params_minimap2 = ''
} else {
    params_minimap2 = params.params_minimap2
}

if (params.params_samtools_view == true) {
    params_samtools_view = ''
} else {
    params_samtools_view = params.params_samtools_view
}

// Step 3. Print inputs and help
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

    usage: nexus run --nf-workflow long_read_rna_variant_calling_de_souza.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'fastq_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
        --reference_genome_fasta_dict_file  :   Reference genome DICT file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.dict).
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
        --deepvariant_input_path            :   DeepVariant input path (default: '/datastore/').
        --deepvariant_output_path           :   DeepVariant output path (default: '/datastore/').
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genome_fasta_fai_file     :   ${params.reference_genome_fasta_fai_file}
        reference_genome_fasta_dict_file    :   ${params.reference_genome_fasta_dict_file}
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
        delete_work_dir                     :   ${params.delete_work_dir}
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
workflow LONG_READ_DNA_VARIANT_CALLING_DE_SOUZA {
    take:
        input_fastq_files_ch             // channel: [val(sample_id), path(fastq_file)]
        reference_genome_fasta_file
        reference_genome_fasta_fai_file
        reference_genome_fasta_dict_file
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
        runMinimap2(
            input_fastq_files_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            params_minimap2,
            platform_tag,
            platform_unit_tag,
            library_tag
        )
        runSamtoolsSamToBam(
            runMinimap2.out.f,
        )
        runSamtoolsFilter(
            runSamtoolsSamToBam.out.f,
            params_samtools_view
        )
        runSamtoolsIndex1(
            runSamtoolsFilter.out.f
        )
        runSamtoolsSort(
            runSamtoolsIndex1.out.f
        )
        run_samtools_sort_output_ch = runSamtoolsSort.out.f
        runGatk4SplitNCigarReads(
            run_samtools_sort_output_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            reference_genome_fasta_dict_file
        )
        run_gatk4_splitncigarreads_output_ch = runGatk4SplitNCigarReads.out.f
        runSamtoolsIndex2(
            run_gatk4_splitncigarreads_output_ch
        )
        run_samtools_index_output_ch = runSamtoolsIndex2.out.f
        run_flag_correction_input_ch = run_samtools_sort_output_ch.join(run_samtools_index_output_ch)
            .map { sample_id, bam_file, bam_bai_file, splitncigarreads_bam_file, splitncigarreads_bam_bai_file ->
                [sample_id, bam_file, bam_bai_file, splitncigarreads_bam_file, splitncigarreads_bam_bai_file]
            }
        runFlagCorrection(
            run_flag_correction_input_ch
        )
        runSamtoolsIndex3(
            runFlagCorrection.out.f
        )
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
                reference_genome_fasta_file,
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
                reference_genome_fasta_file,
                deepvariant_model_type,
                deepvariant_bin_version,
                deepvariant_bin_path,
                deepvariant_input_path,
                deepvariant_output_path,
                output_dir
            )
        }
}

workflow {
    LONG_READ_DNA_VARIANT_CALLING_DE_SOUZA(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_fai_file,
        params.reference_genome_fasta_dict_file,
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

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
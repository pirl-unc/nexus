#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runHLAProfilerPredict } from '../../modules/hlaprofiler'

// Step 2. Input parameters
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.params_hlaprofiler = '-allele_refinement all -if'
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         =================================================================================
         Profile HLA alleles using paired-end RNA sequencing FASTQ files using HLAProfiler
         =================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Profile HLA alleles using paired-end RNA sequencing FASTQ files using HLAProfiler.

    usage: nexus run --nf-workflow paired-end_read_rna_hla_typing_hlaprofiler.nf [required] [optional] [--help]

    required arguments:
        -c                              :   Nextflow .config file.
        -w                              :   Nextflow work directory path.
        --samples_tsv_file              :   TSV file with the following columns:
                                            'sample_id', 'fastq_file_1', 'fastq_file_2'.
        --output_dir                    :   Directory to which output files will be copied.

    optional arguments:
        --params_hlaprofiler            :   HLAProfiler 'predict' parameters (default: '"-allele_refinement all -if"').
                                            Note that the parameters need to be wrapped in quotes.
        --delete_work_dir               :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        params_hlaprofiler              :   ${params.params_hlaprofiler}
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

// Step 5. Workflow
workflow PAIRED_END_RNA_READ_HLA_TYPING_HLAPROFILER {
    take:
        input_fastq_files_ch
        params_hlaprofiler
        output_dir
    main:
        runHLAProfilerPredict(
            input_fastq_files_ch,
            params_hlaprofiler,
            output_dir
        )
    emit:
        runHLAProfilerPredict.out.f
}

workflow {
    PAIRED_END_RNA_READ_HLA_TYPING_HLAPROFILER(
        input_fastq_files_ch,
        params.params_hlaprofiler,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
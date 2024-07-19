#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runStringTie2 } from '../../modules/stringtie2'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.gtf_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf.gz'
params.params_stringtie2 = '-L'
params.delete_work_dir = false

if (params.params_stringtie2 == true) {
    params_stringtie2 = ''
} else {
    params_stringtie2 = params.params_stringtie2
}

// Step 3. Print inputs and help
log.info """\
         =================================================
         Assemble long-read RNA BAM files using StringTie2
         =================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Assemble transcripts using StringTie2.

    usage: nexus run --nf-workflow transcriptome_assembly_stringtie2.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --gtf_file                          :   GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf.gz).
        --params_stringtie2                 :   StringTie2 parameters (default: '"-L"').
                                                Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        gtf_file                            :   ${params.gtf_file}
        params_stringtie2                   :   ${params_stringtie2}
        delete_work_dir                     :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow TRANSCRIPTOME_ASSEMBLY_STRINGTIE2 {
    take:
        input_bam_files_ch            // channel: [val(sample_id), path(bam_file)]
        gtf_file
        params_stringtie2
        output_dir
    main:
        runStringTie2(
            input_bam_files_ch,
            gtf_file,
            params_stringtie2,
            output_dir
        )
    emit:
        runStringTie2.out.f
}

workflow {
    TRANSCRIPTOME_ASSEMBLY_STRINGTIE2(
        input_bam_files_ch,
        params.gtf_file,
        params_stringtie2,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
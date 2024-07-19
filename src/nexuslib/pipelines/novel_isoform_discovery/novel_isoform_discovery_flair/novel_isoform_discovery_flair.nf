//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runFlairAlign } from '../../modules/flair'
include { runFlairCorrect } from '../../modules/flair'
include { runFlairCollapse } from '../../modules/flair'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.gtf_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf'
params.params_flair_align = ''
params.params_flair_correct = ''
params.params_flair_collapse = ''
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         ===================================
         Discover novel isoforms using Flair
         ===================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run flair align.
        2. Run flair correct.
        3. Run flair collapse.

    usage: nexus run --nf-workflow novel_isoform_discovery_flair.nf [required] [optional] [--help]

    required arguments:
        -c                              :   Nextflow .config file.
        -w                              :   Nextflow work directory path.
        --samples_tsv_file              :   TSV file with the following columns:
                                            'sample_id', 'fastq_file'.
        --output_dir                    :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --gtf_file                      :   Reference transcriptome GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf).
        --params_flair_align            :   flair 'align' parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
        --params_flair_correct          :   flair 'correct' parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
        --params_flair_collapse         :   flair 'collapse' parameters (default: '" "').
                                            Note that the parameters need to be wrapped in quotes.
        --delete_work_dir               :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        reference_genome_fasta_file     :   ${params.reference_genome_fasta_file}
        gtf_file                        :   ${params.gtf_file}
        params_flair_align              :   ${params.params_flair_align}
        params_flair_correct            :   ${params.params_flair_correct}
        params_flair_collapse           :   ${params.params_flair_collapse}
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
workflow NOVEL_ISOFORM_DISCOVERY_FLAIR {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        reference_genome_fasta_file
        gtf_file
        params_flair_align
        params_flair_correct
        params_flair_collapse
        output_dir
    main:
        run_flair_align_input_ch = input_fastq_files_ch
        runFlairAlign(
            run_flair_align_input_ch,
            reference_genome_fasta_file,
            params_flair_align,
            output_dir
        )
        runFlairCorrect(
            runFlairAlign.out.f,
            reference_genome_fasta_file,
            gtf_file,
            params_flair_correct,
            output_dir
        )
        runFlairCorrect.out.f.set{ run_flair_correct_output_ch }
        run_flair_correct_output_ch
           .join(input_fastq_files_ch)
           .set{ run_flair_collapse_input_ch }
        runFlairCollapse(
            run_flair_collapse_input_ch,
            reference_genome_fasta_file,
            gtf_file,
            params_flair_collapse,
            output_dir
        )
    emit:
        runFlairCollapse.out.f
}

workflow {
    NOVEL_ISOFORM_DISCOVERY_FLAIR(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.gtf_file,
        params.params_flair_align,
        params.params_flair_correct,
        params.params_flair_collapse,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
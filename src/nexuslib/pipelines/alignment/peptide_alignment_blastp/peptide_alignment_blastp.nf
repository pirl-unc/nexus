#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runBlastp } from '../../modules/blastp'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_proteome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.pc_translations.fa'
params.params_blastp = ''
params.delete_work_dir = false

if (params.params_blastp == true) {
    params_blastp = ''
} else {
    params_blastp = params.params_blastp
}

// Step 3. Print inputs and help
log.info """\
         ======================================
         Align peptide fasta files using Blastp
         ======================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Align peptide sequences (fasta.gz files) to a reference proteome using Blastp.

    usage: nexus run --nf-workflow peptide_alignment_blastp.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'fasta_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_proteome_fasta_file     :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.pc_translations.fa).
        --params_blastp                     :   Blastp parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_proteome_fasta_file       :   ${params.reference_proteome_fasta_file}
        params_blastp                       :   ${params_blastp}
        delete_work_dir                     :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.fasta_file}") }
    .set { input_fasta_files_ch }

// Step 5. Workflow
workflow PEPTIDE_ALIGNMENT_BLASTP {
    take:
        input_fasta_files_ch            // channel: [val(sample_id), path(fasta_file)]
        reference_proteome_fasta_file
        params_blastp
        output_dir
    main:
        runBlastp(
            input_fasta_files_ch,
            reference_proteome_fasta_file,
            params_blastp,
            output_dir
        )
    emit:
        runBlastp.out.f
}

workflow {
    PEPTIDE_ALIGNMENT_BLASTP(
        input_fasta_files_ch,
        params.reference_proteome_fasta_file,
        params_blastp,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

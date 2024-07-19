#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runSvimAlignmentMode } from '../../modules/svim'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.params_svim = '--min_mapq 20 --min_sv_size 30 --insertion_sequences --read_names --zmws'
params.delete_work_dir = false

if (params.params_svim == true) {
    params_svim = ''
} else {
    params_svim = params.params_svim
}

// Step 3. Print inputs and help
log.info """\
         =============================================================================
         Identify structural variants in long-read DNA sequencing BAM files using SVIM
         =============================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run SVIM.

    usage: nexus run --nf-workflow long_read_dna_variant_calling_svim.nf [required] [optional] [--help]

    required arguments:
        -c                              :   Nextflow .config file.
        -w                              :   Nextflow work directory path.
        --samples_tsv_file              :   TSV file with the following columns:
                                            'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                    :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --params_svim                   :   SVIM parameters (default: '"--min_mapq 20 --min_sv_size 30 --insertion_sequences --read_names --zmws"').
                                            Note that the parameters need to be wrapped in quotes.
        --delete_work_dir               :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        reference_genome_fasta_file     :   ${params.reference_genome_fasta_file}
        params_svim                     :   ${params_svim}
        delete_work_dir                 :   ${params.delete_work_dir}
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
workflow LONG_READ_DNA_VARIANT_CALLING_SVIM {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        params_svim
        output_dir

    main:
        runSvimAlignmentMode(
            input_bam_files_ch,
            reference_genome_fasta_file,
            params_svim,
            output_dir
        )
}

workflow {
    LONG_READ_DNA_VARIANT_CALLING_SVIM(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params_svim,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

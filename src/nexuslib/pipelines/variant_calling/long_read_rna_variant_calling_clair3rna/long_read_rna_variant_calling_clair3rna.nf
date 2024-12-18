#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runClair3RNA } from '../../modules/clair3rna'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.params_clair3rna = '--platform hifi_sequel2_minimap2'
params.delete_work_dir = false

if (params.params_clair3rna == true) {
    params_clair3rna = ''
} else {
    params_clair3rna = params.params_clair3rna
}

// Step 3. Print inputs and help
log.info """\
         =========================================================================================================
         Identify RNA variants in long-read RNA sequencing BAM files using Clair3-RNA (Zheng et al., bioRxiv 2024)
         =========================================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run Clair3-RNA.

    usage: nexus run --nf-workflow long_read_rna_variant_calling_clair3rna.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
        --params_clair3rna                  :   Clair3-RNA parameters (default: '"--platform hifi_sequel2_minimap2"').
                                                Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genome_fasta_fai_file     :   ${params.reference_genome_fasta_fai_file}
        params_clair3rna                    :   ${params_clair3rna}
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
workflow LONG_READ_DNA_VARIANT_CALLING_CLAIR3RNA {
    take:
        input_bam_files_ch              // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        reference_genome_fasta_fai_file
        params_clair3rna
        output_dir

    main:
        runClair3RNA(
            input_bam_files_ch,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            params_clair3rna,
            output_dir
        )
}

workflow {
    LONG_READ_DNA_VARIANT_CALLING_CLAIR3RNA(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_fai_file,
        params_clair3rna,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

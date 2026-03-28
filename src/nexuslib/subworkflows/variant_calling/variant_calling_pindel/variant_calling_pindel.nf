#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }                  from '../../../tools/samtools'
include { runPindel }                              from '../../../tools/pindel'
include { decompressFile as decompressFasta }      from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help = ''

// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.reference_genome_fasta_file              = ''

// Optional arguments
params.params_pindel                            = ''

if (params.params_pindel == true) {
    params_pindel = ''
} else {
    params_pindel = params.params_pindel
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ===================================================================================
         Identify germline variants in paired-end read DNA sequencing BAM files using Pindel
         ===================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run Pindel.

    usage: nexus run --nf-workflow variant_calling_pindel.nf [required] [optional] [--help]

    required arguments:
        --samples_tsv_file                      :   TSV file with the following columns:
                                                    'sample_id',
                                                    'bam_file',
                                                    'bam_bai_file'.
        --output_dir                            :   Directory to which output files will be symlinked.
        --reference_genome_fasta_file           :   Reference genome FASTA file.

    optional arguments:
        --params_pindel                         :   Pindel parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                        :   ${params.samples_tsv_file}
        output_dir                              :   ${params.output_dir}
        reference_genome_fasta_file             :   ${params.reference_genome_fasta_file}
        params_pindel                           :   ${params_pindel}
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
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_PINDEL {
    take:
        input_bam_files_ch          // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        output_dir
        reference_genome_fasta_file
        params_pindel

    main:
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidxFasta(decompressFasta.out.f)
        fasta_file      = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file  = runSamtoolsFaidxFasta.out.fai_file

        runPindel(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            params_pindel,
            output_dir
        )

    emit:
        runPindel.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_PINDEL(
        input_bam_files_ch,
        params.output_dir,
        params.reference_genome_fasta_file,
        params_pindel
    )
}

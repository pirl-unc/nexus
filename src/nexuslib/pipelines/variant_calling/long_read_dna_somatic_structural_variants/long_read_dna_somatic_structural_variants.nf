#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runSavana } from '../../../../modules/savana'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.reference_genome_fasta_file = ''
params.savana = ''
params.savana_params = ''
// Optional arguments
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
        =============================================================================================================
        Identify Somatic Structural DNA Variants Using Pacific Biosciences HiFi CCS Whole-genome Sequencing BAM Files
        =============================================================================================================
        PURPOSE:
            This workflow is intended for long-read DNA sequencing BAM files.

        WORKFLOW:
            1. Run Savana.
         """.stripIndent()

if (params.help) {
    log.info"""\
        usage: nextflow run pacbio_dna_somatic_structural_variants.nf [required] [optional] [--help]

        required arguments:
            --samples_tsv_file              :   TSV file with the following columns:
                                                'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'
            --output_dir                    :   Directory to which output files will be copied.
            --reference_genome_fasta_file   :   Reference genome FASTA file.
            --savana                        :   Savana path.
            --savana_params                 :   Savana parameters (e.g. "").

        optional arguments:
            --delete_work_dir               :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
            samples_tsv_file                :   ${params.samples_tsv_file}
            output_dir                      :   ${params.output_dir}
            reference_genome_fasta_file     :   ${params.reference_genome_fasta_file}
            savana                          :   ${params.savana}
            savana_params                   :   ${params.savana_params}
            delete_work_dir                 :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.tumor_bam_file}",
        "${row.tumor_bam_bai_file}",
        "${row.normal_bam_file}",
        "${row.normal_bam_bai_file}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow PACBIO_DNA_SOMATIC_STRUCTURAL_VARIANTS {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)]
        reference_genome_fasta_file
        savana
        savana_params
        output_dir

    main:
        run_savana_input_ch = input_bam_files_ch

        // Savana
        runSavana(
            run_savana_input_ch,
            reference_genome_fasta_file,
            savana,
            savana_params,
            output_dir
        )
        runSavana.out.f.set{ run_savana_output_ch }
}

workflow {
    PACBIO_DNA_SOMATIC_STRUCTURAL_VARIANTS(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.savana,
        params.savana_params,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

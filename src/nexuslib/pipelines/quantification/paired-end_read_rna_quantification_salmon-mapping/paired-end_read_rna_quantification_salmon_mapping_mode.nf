#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 * Last updated date: Oct 3, 2023
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runSalmonPairedEndMappingMode } from '../../../modules/salmon'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.salmon = ''
params.salmon_index_params = ''
params.salmon_quant_params = ''
// Optional arguments
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         ===================================================================
         Quantify RNA Abundance with FASTQ Files Using Salmon (Mapping Mode)
         ===================================================================
         PURPOSE:
            This workflow is intended for paired-end RNA sequencing FASTQ files.

         WORKFLOW:
            1. Run Salmon (mapping mode).
         """.stripIndent()

if (params.help) {
    log.info"""\
        usage: nextflow run paired_end_read_rna_quantification_salmon_mapping_mode.nf [required] [optional] [--help]

        required arguments:
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'fastq_file_1', 'fastq_file_2', 'fasta_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --salmon                            :   salmon path.
            --salmon_index_params               :   salmon 'index' parameters (e.g. '--gencode').
            --salmon_quant_params               :   salmon 'quant' parameters
                                                    (e.g. '--libType <library_type> --geneMap <gtf_file> --seqBias --gcBias --posBias').
            --salmon_index                      :   salmon index path.

        optional arguments:
            --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
            samples_tsv_file                    :   ${params.samples_tsv_file}
            output_dir                          :   ${params.output_dir}
            salmon                              :   ${params.salmon}
            salmon_index_params                 :   ${params.salmon_index_params}
            salmon_quant_params                 :   ${params.salmon_quant_params}
            delete_work_dir                     :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.fastq_file_1}",
        "${row.fastq_file_2}",
        "${row.fasta_file}") }
    .set { input_fastq_files_ch }

// Step 5. Workflow
workflow PAIRED_END_RNA_READ_QUANTIFICATION_SALMON_MAPPING_MODE {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file_1), path(fastq_file_2), path(fasta_file)]
        salmon
        salmon_index_params
        salmon_quant_params
        output_dir
    main:
        run_salmon_input_ch = input_fastq_files_ch
        runSalmonPairedEndMappingMode(
            run_salmon_input_ch,
            salmon,
            salmon_index_params,
            salmon_quant_params,
            output_dir
        )
    emit:
        runSalmonPairedEndMappingMode.out.f
}

workflow {
    PAIRED_END_RNA_READ_QUANTIFICATION_SALMON_MAPPING_MODE(
        input_fastq_files_ch,
        params.salmon,
        params.salmon_index_params,
        params.salmon_quant_params,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
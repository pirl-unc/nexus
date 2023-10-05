#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 * Last updated date: Oct 3, 2023
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runDelly2TumorNormal } from '../../../modules/delly2'
include { runLumpyExpressTumorNormal } from '../../../modules/lumpy'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.reference_genome_fasta_file = ''
params.delly2 = ''
params.delly2_call_params = ''
params.bcftools = ''
params.python2 = ''
params.lumpy_express = ''
params.lumpy_extract_split_reads_script_file = ''
params.lumpy_config_file = ''
params.samtools = ''
// Optional arguments

// Step 3. Print inputs and help
log.info """\
         =========================================================================================
         Identify Somatic Structural Variants using (Illumina) Paired-end DNA Sequencing BAM Files
         =========================================================================================
         PURPOSE:
            This workflow is intended for pairs of tumor and matched normal samples.

         IMPORTANT NOTES:
            - Please run this script in a python2 environment.

         WORKFLOW:
            1. Run DELLY2 (tumor and normal mode).
            2. Run LUMPY Express (tumor and normal mode).
         """.stripIndent()

if (params.help) {
    log.info"""\
        usage: nextflow run paired_end_dna_somatic_structural_variants.nf [required] [optional] [--help]

        required arguments:
            --samples_tsv_file                          :   TSV file with the following columns:
                                                            'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file', 'tumor_sample_id', 'normal_sample_id'
            --output_dir                                :   Directory to which output files will be symlinked.
            --reference_genome_fasta_file               :   Reference genome FASTA file.
            --delly2                                    :   DELLY2 path.
            --delly2_call_params                        :   DELLY2 'call' parameters (e.g. "").
            --bcftools                                  :   bcftools path.
            --python2                                   :   python2 path.
            --lumpy_express                             :   LUMPY Express path.
            --lumpy_extract_split_reads_script_file     :   LUMPY 'extractSplitReads_BwaMem' file.
            --lumpy_config_file                         :   LUMPY config file.
            --samtools                                  :   samtools path.

        optional arguments:
    """.stripIndent()
    exit 0
} else {
    log.info"""\
            samples_tsv_file                            :   ${params.samples_tsv_file}
            output_dir                                  :   ${params.output_dir}
            reference_genome_fasta_file                 :   ${params.reference_genome_fasta_file}
            delly2                                      :   ${params.delly2}
            delly2_call_params                          :   ${params.delly2_call_params}
            bcftools                                    :   ${params.bcftools}
            python2                                     :   ${params.python2}
            lumpy_express                               :   ${params.lumpy_express}
            lumpy_extract_split_reads_script_file       :   ${params.lumpy_extract_split_reads_script_file}
            lumpy_config_file                           :   ${params.lumpy_config_file}
            samtools                                    :   ${params.samtools}
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
        "${row.normal_bam_bai_file}",
        "${row.tumor_sample_id}",
        "${row.normal_sample_id}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow PAIRED_END_DNA_SOMATIC_STRUCTURAL_VARIANTS {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(tumor_sample_id), val(normal_sample_id)]
        output_dir
        reference_genome_fasta_file
        delly2
        delly2_call_params
        bcftools
        python2
        lumpy_express
        lumpy_extract_split_reads_script_file
        lumpy_config_file
        samtools

    main:
        input_bam_files_ch
            .map{ [it[0], it[1], it[2], it[3], it[4], it[5], it[6]] }
            .set{ run_delly2_input_ch }
        input_bam_files_ch
            .map{ [it[0], it[1], it[2], it[3], it[4], it[5], it[6]] }
            .set{ run_lumpy_input_ch }

        // DELLY2
        runDelly2TumorNormal(
            run_delly2_input_ch,
            delly2,
            delly2_call_params,
            bcftools,
            reference_genome_fasta_file,
            output_dir
        )

        // LUMPY Express
        runLumpyExpressTumorNormal(
            run_lumpy_input_ch,
            python2,
            lumpy_express,
            lumpy_extract_split_reads_script_file,
            lumpy_config_file,
            samtools,
            output_dir
        )
}

workflow {
    PAIRED_END_DNA_SOMATIC_STRUCTURAL_VARIANTS(
        input_bam_files_ch,
        params.output_dir,
        params.reference_genome_fasta_file,
        params.delly2,
        params.delly2_call_params,
        params.bcftools,
        params.python2,
        params.lumpy_express,
        params.lumpy_extract_split_reads_script_file,
        params.lumpy_config_file,
        params.samtools
    )
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 * Last updated date: Oct 4, 2023
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runSniffles2 } from '../../../../modules/sniffles2'
include { runPbsv } from '../../../../modules/pbsv'
include { runSvimAlignmentMode } from '../../../../modules/svim'
include { runCuteSV } from '../../../../modules/cutesv'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.reference_genome_fasta_file = ''
params.sniffles2 = ''
params.sniffles2_params = ''
params.pbsv = ''
params.pbsv_discover_params = ''
params.pbsv_call_params = ''
params.svim = ''
params.svim_params = ''
params.cutesv_params = ''
// Optional arguments
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
        =====================================================================================================
        Identify Structural DNA Variants Using Pacific Biosciences HiFi CCS Whole-genome Sequencing BAM Files
        =====================================================================================================
        PURPOSE:
            This workflow is intended for long-read DNA sequencing BAM files.

        WORKFLOW:
            1. Run Sniffles2, PBSV, SVIM, and cuteSV.
         """.stripIndent()

if (params.help) {
    log.info"""\
        usage: nextflow run pacbio_dna_structural_variants.nf [required] [optional] [--help]

        required arguments:
            --samples_tsv_file              :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
            --output_dir                    :   Directory to which output files will be copied.
            --reference_genome_fasta_file   :   Reference genome FASTA file.
            --sniffles2                     :   Sniffles2 path.
            --sniffles2_params              :   Sniffles2 parameters (e.g. "").
            --pbsv                          :   pbsv path.
            --pbsv_discover_params          :   pbsv 'discover' parameters (e.g. "").
            --pbsv_call_params              :   pbsv 'call' parameters (e.g. "").
            --svim                          :   SVIM path.
            --svim_params                   :   SVIM parameters (e.g. "").
            --cutesv                        :   cuteSV path.
            --cutesv_params                 :   cuteSV parameters (e.g. "").

        optional arguments:
            --delete_work_dir               :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
            samples_tsv_file                :   ${params.samples_tsv_file}
            output_dir                      :   ${params.output_dir}
            reference_genome_fasta_file     :   ${params.reference_genome_fasta_file}
            sniffles2                       :   ${params.sniffles2}
            sniffles2_params                :   ${params.sniffles2_params}
            pbsv                            :   ${params.pbsv}
            pbsv_discover_params            :   ${params.pbsv_discover_params}
            pbsv_call_params                :   ${params.pbsv_call_params}
            svim                            :   ${params.svim}
            svim_params                     :   ${params.svim_params}
            cutesv                          :   ${params.cutesv}
            cutesv_params                   :   ${params.cutesv_params}
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
workflow PACBIO_DNA_STRUCTURAL_VARIANTS {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        sniffles2
        sniffles2_params
        pbsv
        pbsv_discover_params
        pbsv_call_params
        svim
        svim_params
        cutesv
        cutesv_params
        output_dir

    main:
        run_sniffles2_input_ch = input_bam_files_ch
        run_pbsv_input_ch = input_bam_files_ch
        run_svim_input_ch = input_bam_files_ch
        run_cutesv_input_ch = input_bam_files_ch

        // Sniffles2
        runSniffles2(
            run_sniffles2_input_ch,
            reference_genome_fasta_file,
            sniffles2,
            sniffles2_params,
            output_dir
        )
        runSniffles2.out.f.set{ run_sniffles2_output_ch }

        // pbsv
        runPbsv(
            run_pbsv_input_ch,
            reference_genome_fasta_file,
            pbsv,
            pbsv_discover_params,
            pbsv_call_params,
            output_dir
        )
        runPbsv.out.f.set{ run_pbsv_output_ch }

        // SVIM
        runSvimAlignmentMode(
            run_svim_input_ch,
            reference_genome_fasta_file,
            svim,
            svim_params,
            output_dir
        )
        runSvimAlignmentMode.out.f.set{ run_svim_output_ch }

        // cuteSV
        runCuteSV(
            run_cutesv_input_ch,
            reference_genome_fasta_file,
            cutesv,
            cutesv_params,
            output_dir
        )
        runCuteSV.out.f.set{ run_cutesv_output_ch }
}

workflow {
    PACBIO_DNA_STRUCTURAL_VARIANTS(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.sniffles2,
        params.sniffles2_params,
        params.pbsv,
        params.pbsv_discover_params,
        params.pbsv_call_params,
        params.svim,
        params.svim_params,
        params.cutesv,
        params.cutesv_params,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

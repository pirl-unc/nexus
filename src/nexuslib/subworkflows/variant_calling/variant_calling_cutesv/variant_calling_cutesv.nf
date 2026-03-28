#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }   from '../../../tools/samtools'
include { runCuteSV }               from '../../../tools/cutesv'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir                   = ''
params.reference_genome_fasta_file  = ''

// Optional arguments
params.params_cutesv                = '--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_support 3 --min_mapq 20 --min_size 30 --max_size -1 --report_readid --genotype'

if (params.params_cutesv == true) {
    params_cutesv = ''
} else {
    params_cutesv = params.params_cutesv
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ===============================================================================
         Identify structural variants in long-read DNA sequencing BAM files using CuteSV
         ===============================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run CuteSV.

    usage: nexus run --nf-workflow variant_calling_cutesv.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.

    optional arguments:
        --params_cutesv                     :   CuteSV parameters (default:
                                                '"--max_cluster_bias_INS 1000
                                                  --diff_ratio_merging_INS 0.9
                                                  --max_cluster_bias_DEL 1000
                                                  --diff_ratio_merging_DEL 0.5
                                                  --min_support 3
                                                  --min_mapq 20
                                                  --min_size 30
                                                  --max_size -1
                                                  --report_readid
                                                  --genotype"').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_cutesv                       :   ${params_cutesv}
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
workflow VARIANT_CALLING_CUTESV {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        params_cutesv
        output_dir

    main:
        runSamtoolsFaidxFasta(reference_genome_fasta_file)
        fasta_file      = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file  = runSamtoolsFaidxFasta.out.fai_file

        runCuteSV(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            params_cutesv,
            output_dir
        )

    emit:
        runCuteSV.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_CUTESV(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params_cutesv,
        params.output_dir
    )
}

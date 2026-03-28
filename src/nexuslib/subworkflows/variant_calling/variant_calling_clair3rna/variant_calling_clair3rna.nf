#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }    from '../../../tools/samtools'
include { runClair3RNA }        from '../../../tools/clair3rna'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''

// Optional arguments
params.params_clair3rna                 = '--platform hifi_sequel2_minimap2'

if (params.params_clair3rna == true) {
    params_clair3rna = ''
} else {
    params_clair3rna = params.params_clair3rna
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         =========================================================================================================
         Identify RNA variants in long-read RNA sequencing BAM files using Clair3-RNA (Zheng et al., bioRxiv 2024)
         =========================================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run Clair3-RNA.

    usage: nexus run --nf-workflow variant_calling_clair3rna.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.

    optional arguments:
        --params_clair3rna                  :   Clair3-RNA parameters (default: '"--platform hifi_sequel2_minimap2"').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_clair3rna                    :   ${params_clair3rna}
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
workflow VARIANT_CALLING_CLAIR3RNA {
    take:
        input_bam_files_ch              // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        params_clair3rna
        output_dir

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file      = runSamtoolsFaidx.out.fasta
        fasta_fai_file  = runSamtoolsFaidx.out.fai_file

        runClair3RNA(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            params_clair3rna,
            output_dir
        )

    emit:
        runClair3RNA.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_CLAIR3RNA(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params_clair3rna,
        params.output_dir
    )
}


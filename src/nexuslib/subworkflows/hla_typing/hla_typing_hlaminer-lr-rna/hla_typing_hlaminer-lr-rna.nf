#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runHLAminerLongReadRNA }   from '../../../tools/hlaminer'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_minimap2      = '-ax map-hifi --secondary=no'
params.params_hlaminer      = '-e 1 -s 500 -q 1 -i 1'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow HLA_TYPING_HLAMINER_LR_RNA {
    take:
        input_fastq_files_ch
        params_minimap2
        params_hlaminer
        output_dir

    main:
        runHLAminerLongReadRNA(
            input_fastq_files_ch,
            params_minimap2,
            params_hlaminer,
            output_dir
        )
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ==================================================================
             Profile HLA alleles from long-read RNA-seq FASTQ files using
             HLAminer (minimap2 alignment against HLAminer's HLA CDS reference)
             ==================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Align long-read RNA-seq FASTQ files against HLAminer's HLA CDS
               reference using minimap2.
            2. Pipe the alignment SAM straight into HLAminer.pl to predict HLA
               class I and II alleles.

        usage: nexus run --nf-workflow hla_typing_hlaminer-lr-rna.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_minimap2       :   minimap2 extra CLI parameters
                                        (default: '"-ax map-hifi --secondary=no"').
                                        Use -ax map-hifi for PacBio HiFi reads,
                                        or -ax map-ont for Oxford Nanopore reads.
                                        Note that the parameters need to be wrapped in quotes.
            --params_hlaminer       :   HLAminer.pl extra CLI parameters
                                        (default: '"-e 1 -s 500 -q 1 -i 1"').
                                        Long reads need relaxed thresholds; the
                                        default mirrors HLAminer's official
                                        long-read RNA-seq demo (HPRArnaseq_ONT*).
                                        Common flags include
                                          -i <minimum percent identity>
                                          -s <minimum alignment score>
                                          -q <minimum log10 expect value>
                                          -e <single-end reads (1=yes/0=no)>
                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_minimap2 = (params.params_minimap2 == true) ? '' : params.params_minimap2
    def params_hlaminer = (params.params_hlaminer == true) ? '' : params.params_hlaminer

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_minimap2         :   ${params_minimap2}
        params_hlaminer         :   ${params_hlaminer}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_fastq_files_ch }

    HLA_TYPING_HLAMINER_LR_RNA(
        input_fastq_files_ch,
        params_minimap2,
        params_hlaminer,
        params.output_dir
    )
}

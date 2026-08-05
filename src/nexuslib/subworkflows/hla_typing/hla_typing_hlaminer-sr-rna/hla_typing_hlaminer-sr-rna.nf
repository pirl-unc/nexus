#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runHLAminerShortReadRNA }   from '../../../tools/hlaminer'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_bwa_aln       = '-e 0 -o 0'
params.params_bwa_sampe     = '-o 1000'
params.params_hlaminer      = '-s 500'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow HLA_TYPING_HLAMINER_SR_RNA {
    take:
        input_fastq_files_ch
        params_bwa_aln
        params_bwa_sampe
        params_hlaminer
        output_dir

    main:
        runHLAminerShortReadRNA(
            input_fastq_files_ch,
            params_bwa_aln,
            params_bwa_sampe,
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
             Profile HLA alleles from short-read RNA FASTQ files using HLAminer
             (bwa alignment against HLAminer's coding HLA reference)
             ==================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Align paired short-read RNA-seq FASTQ files against HLAminer's
               coding HLA reference (HLA-I_II_CDS.fasta) using bwa aln/sampe.
            2. Run HLAminer.pl in file mode on the paired SAM to predict HLA
               class I and II alleles.

        usage: nexus run --nf-workflow hla_typing_hlaminer-sr-rna.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file_1', 'fastq_file_2'.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_bwa_aln        :   bwa aln extra CLI parameters
                                        (default: '"-e 0 -o 0"').
                                        Note that the parameters need to be wrapped in quotes.
            --params_bwa_sampe      :   bwa sampe extra CLI parameters
                                        (default: '"-o 1000"').
                                        Note that the parameters need to be wrapped in quotes.
            --params_hlaminer       :   HLAminer.pl extra CLI parameters
                                        (default: '"-s 500"').
                                        Mirrors HLAminer's official short-read demo
                                        (HPRArnaseq_classI-II.sh); short reads keep the
                                        strict default identity/expect thresholds.
                                        Common flags include
                                          -i <minimum percent identity>
                                          -s <minimum alignment score>
                                          -q <minimum log10 expect value>
                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_bwa_aln = (params.params_bwa_aln == true) ? '' : params.params_bwa_aln
    def params_bwa_sampe = (params.params_bwa_sampe == true) ? '' : params.params_bwa_sampe
    def params_hlaminer = (params.params_hlaminer == true) ? '' : params.params_hlaminer

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_bwa_aln          :   ${params_bwa_aln}
        params_bwa_sampe        :   ${params_bwa_sampe}
        params_hlaminer         :   ${params_hlaminer}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file_1}",
            "${row.fastq_file_2}") }
        .set { input_fastq_files_ch }

    HLA_TYPING_HLAMINER_SR_RNA(
        input_fastq_files_ch,
        params_bwa_aln,
        params_bwa_sampe,
        params_hlaminer,
        params.output_dir
    )
}

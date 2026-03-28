//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runEspresso }                        from '../../../tools/espresso'
include { decompressFile as decompressFasta }  from '../../../tools/utils'
include { decompressFile as decompressGtf }    from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir                   = ''
params.reference_genome_fasta_file  = ''
params.reference_genes_gtf_file     = ''

// Optional arguments
params.params_espresso_s            = ''
params.params_espresso_c            = ''
params.params_espresso_q            = ''

if (params.params_espresso_s == true) {
    params_espresso_s = ''
} else {
    params_espresso_s = params.params_espresso_s
}

if (params.params_espresso_c == true) {
    params_espresso_c = ''
} else {
    params_espresso_c = params.params_espresso_c
}

if (params.params_espresso_q == true) {
    params_espresso_q = ''
} else {
    params_espresso_q = params.params_espresso_q
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ====================================
         Characterize isoforms using ESPRESSO
         ====================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run ESPRESSO.

    usage: nexus run --nf-workflow isoform_characterization_espresso.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file (uncompressed .fasta file).
        --reference_genes_gtf_file          :   Reference genes GTF file (uncompressed .gtf file).

    optional arguments:
        --params_espresso_s                 :   ESPRESSO_S.pl parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
        --params_espresso_c                 :   ESPRESSO_C.pl parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
        --params_espresso_q                 :   ESPRESSO_Q.pl parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_espresso_s                   :   ${params_espresso_s}
        params_espresso_c                   :   ${params_espresso_c}
        params_espresso_q                   :   ${params_espresso_q}
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
workflow ISOFORM_CHARACTERIZATION_ESPRESSO {
    take:
        input_bam_files_ch
        reference_genome_fasta_file
        reference_genes_gtf_file
        params_espresso_s
        params_espresso_c
        params_espresso_q
        output_dir

    main:
        decompressFasta(reference_genome_fasta_file)
        decompressGtf(reference_genes_gtf_file)

        runEspresso(
            input_bam_files_ch,
            decompressFasta.out.f,
            decompressGtf.out.f,
            params_espresso_s,
            params_espresso_c,
            params_espresso_q,
            output_dir
        )

    emit:
        runEspresso.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ISOFORM_CHARACTERIZATION_ESPRESSO(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params_espresso_s,
        params_espresso_c,
        params_espresso_q,
        params.output_dir
    )
}

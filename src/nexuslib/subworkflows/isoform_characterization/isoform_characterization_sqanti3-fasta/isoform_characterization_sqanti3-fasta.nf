//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSqanti3FastaMode }                    from '../../../tools/sqanti3'
include { decompressFile as decompressFasta }      from '../../../tools/utils'
include { decompressFile as decompressGtf }        from '../../../tools/utils'
include { decompressFile as decompressInput }      from '../../../tools/utils'

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
params.params_sqanti3_qc            = '--aligner_choice minimap2 --report html --force_id_ignore'
params.params_sqanti3_filter        = 'ml'

if (params.params_sqanti3_qc == true) {
    params_sqanti3_qc = ''
} else {
    params_sqanti3_qc = params.params_sqanti3_qc
}

if (params.params_sqanti3_filter == true) {
    params_sqanti3_filter = ''
} else {
    params_sqanti3_filter = params.params_sqanti3_filter
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ================================================
         Characterize isoforms using Sqanti3 (FASTA mode)
         ================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run SQANTI3 (FASTA mode).

    usage: nexus run --nf-workflow isoform_characterization_sqanti3-fasta.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'fasta_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.
        --reference_genes_gtf_file          :   Reference genes GTF file.

    optional arguments:
        --params_sqanti3_qc                 :   sqanti3_qc.py parameters (default: '"--aligner_choice minimap2 --report html --force_id_ignore --isoAnnotLite"').
                                                Note that the parameters need to be wrapped in quotes.
        --params_sqanti3_filter             :   sqanti3_filter.py parameters (default: '"ml"').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_sqanti3_qc                   :   ${params_sqanti3_qc}
        params_sqanti3_filter               :   ${params_sqanti3_filter}
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
        "${row.fasta_file}") }
    .set { input_fasta_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow ISOFORM_CHARACTERIZATION_SQANTI3_FASTA {
    take:
        input_fasta_files_ch
        reference_genome_fasta_file
        reference_genes_gtf_file
        params_sqanti3_qc
        params_sqanti3_filter
        method
        output_dir

    main:
        decompressFasta(reference_genome_fasta_file)
        decompressGtf(reference_genes_gtf_file)

        // Decompress sample FASTA files if needed
        input_fasta_files_ch.map { sid, fa -> fa }.set { input_fasta_raw_ch }
        decompressInput(input_fasta_raw_ch)
        input_fasta_files_ch.map { sid, fa -> sid }
            .merge(decompressInput.out.f)
            .set { input_fasta_decompressed_ch }

        runSqanti3FastaMode(
            input_fasta_decompressed_ch,
            decompressFasta.out.f,
            decompressGtf.out.f,
            params_sqanti3_qc,
            params_sqanti3_filter,
            method,
            output_dir
        )

    emit:
        runSqanti3FastaMode.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ISOFORM_CHARACTERIZATION_SQANTI3_FASTA(
        input_fasta_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params_sqanti3_qc,
        params_sqanti3_filter,
        "sqanti3_fasta",
        params.output_dir
    )
}

//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }    from '../../../tools/samtools'
include { runIsoquant }         from '../../../tools/isoquant'

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
params.params_isoquant              = '--data_type pacbio_ccs --sqanti_output --high_memory --complete_genedb'

if (params.params_isoquant == true) {
    params_isoquant = ''
} else {
    params_isoquant = params.params_isoquant
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ====================================
         Characterize isoforms using Isoquant
         ====================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run isoquant.

    usage: nexus run --nf-workflow isoform_characterization_isoquant.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'fastq_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.
        --reference_genes_gtf_file          :   Reference genes GTF file.

    optional arguments:
        --params_isoquant                   :   isoquant parameters (default: '"--data_type pacbio_ccs --sqanti_output --high_memory --complete_genedb"').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_isoquant                     :   ${params_isoquant}
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
        "${row.fastq_file}") }
    .set { input_fastq_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow ISOFORM_CHARACTERIZATION_ISOQUANT {
    take:
        input_fastq_files_ch
        reference_genome_fasta_file
        reference_genes_gtf_file
        params_isoquant
        output_dir

    main:
        // Step 1. Index reference genome FASTA file
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file          = runSamtoolsFaidx.out.fasta
        fasta_fai_file      = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file      = runSamtoolsFaidx.out.gzi_file

        // Step 2. Run Isoquant
        runIsoquant(
            input_fastq_files_ch,
            fasta_file,
            fasta_fai_file,
            reference_genes_gtf_file,
            params_isoquant,
            output_dir
        )

    emit:
        runIsoquant.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ISOFORM_CHARACTERIZATION_ISOQUANT(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params_isoquant,
        params.output_dir
    )
}

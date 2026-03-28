#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runKallistoIndex }                from '../../../tools/kallisto'
include { runKallistoT2G }                  from '../../../tools/kallisto'
include { runKallistoQuantTccLongReads }    from '../../../tools/kallisto'
include { decompressFile as decompressFasta } from '../../../tools/utils'
include { decompressFile as decompressGtf }   from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                                 = ''

// Required arguments
params.samples_tsv_file                     = ''
params.output_dir                           = ''
params.reference_transcripts_fasta_file     = ''
params.reference_genes_gtf_file             = ''

// Optional arguments
params.params_kallisto_index                = '-k 63'
params.params_kallisto_bus                  = '-x bulk --threshold 0.8'
params.params_bustools_sort                 = ''
params.params_bustools_count                = '--cm -m'
params.params_kallisto_quanttcc             = '-P PacBio --matrix-to-files'

if (params.params_kallisto_index == true) {
    params_kallisto_index = ''
} else {
    params_kallisto_index = params.params_kallisto_index
}

if (params.params_kallisto_bus == true) {
    params_kallisto_bus = ''
} else {
    params_kallisto_bus = params.params_kallisto_bus
}

if (params.params_bustools_sort == true) {
    params_bustools_sort = ''
} else {
    params_bustools_sort = params.params_bustools_sort
}

if (params.params_bustools_count == true) {
    params_bustools_count = ''
} else {
    params_bustools_count = params.params_bustools_count
}

if (params.params_kallisto_quanttcc == true) {
    params_kallisto_quanttcc = ''
} else {
    params_kallisto_quanttcc = params.params_kallisto_quanttcc
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ====================================================
         Quantify RNA in long-read FASTQ files using Kallisto
         ====================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Quantify RNA in long-read FASTQ files using Kallisto (quant-tcc).

    usage: nexus run --nf-workflow quantification_kallisto-lr.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'fastq_file'.
        --reference_transcripts_fasta_file  :   Reference transcripts FASTA file.
        --reference_genes_gtf_file          :   Reference genes GTF file.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --params_kallisto_index             :   Kallisto index parameters (default: '"-k 63"').
        --params_kallisto_bus               :   Kallisto bus parameters (default: '"-x bulk --threshold 0.8"').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
        --params_bustools_sort              :   Bustools sort parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
        --params_bustools_count             :   Bustools count parameters (default: '"--cm -m"').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
        --params_kallisto_quanttcc          :   Kallisto quant-tcc parameters (default: '"-P PacBio --matrix-to-files"').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        reference_transcripts_fasta_file    :   ${params.reference_transcripts_fasta_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_kallisto_index               :   ${params_kallisto_index}
        params_kallisto_bus                 :   ${params_kallisto_bus}
        params_bustools_sort                :   ${params_bustools_sort}
        params_bustools_count               :   ${params_bustools_count}
        params_kallisto_quanttcc            :   ${params_kallisto_quanttcc}
        output_dir                          :   ${params.output_dir}
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
workflow QUANTIFICATION_KALLISTO_LR {
    take:
        input_fastq_files_ch
        reference_transcripts_fasta_file
        reference_genes_gtf_file
        params_kallisto_index
        params_kallisto_bus
        params_bustools_sort
        params_bustools_count
        params_kallisto_quanttcc
        output_dir

    main:
        // Step 1. Decompress reference files if needed
        decompressFasta(reference_transcripts_fasta_file)
        decompressGtf(reference_genes_gtf_file)

        // Step 2. Index reference FASTA file
        runKallistoIndex(
            decompressFasta.out.f,
            params_kallisto_index
        )

        // Step 3. Index reference GTF file
        runKallistoT2G(decompressGtf.out.f)

        // Step 4. Run Kallisto
        runKallistoQuantTccLongReads(
            input_fastq_files_ch,
            runKallistoIndex.out.index_file,
            decompressGtf.out.f,
            runKallistoT2G.out.t2g_file,
            params_kallisto_bus,
            params_bustools_sort,
            params_bustools_count,
            params_kallisto_quanttcc,
            output_dir
        )

    emit:
        runKallistoQuantTccLongReads.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    QUANTIFICATION_KALLISTO_LR(
        input_fastq_files_ch,
        params.reference_transcripts_fasta_file,
        params.reference_genes_gtf_file,
        params_kallisto_index,
        params_kallisto_bus,
        params_bustools_sort,
        params_bustools_count,
        params_kallisto_quanttcc,
        params.output_dir
    )
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runKallistoQuantTccLongReads } from '../../modules/kallisto'

// Step 2. Input parameters
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.kallisto_index_file = '/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/kallisto/v0.51.1/gencode_v41/kallisto_gencode_v41_index_k63.idx'
params.gtf_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf'
params.t2g_file = '/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/kallisto/v0.51.1/gencode_v41/gencode_v41.t2g'
params.params_kallisto_bus = '-x bulk --threshold 0.8'
params.params_bustools_sort = ''
params.params_bustools_count = '--cm -m'
params.params_kallisto_quanttcc = '-P PacBio --matrix-to-files'
// Optional arguments
params.delete_work_dir = false

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

// Step 3. Print inputs and help
log.info """\
         ====================================================
         Quantify RNA in long-read FASTQ files using Kallisto
         ====================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Quantify RNA in long-read FASTQ files using Kallisto (quant-tcc).

    usage: nexus run --nf-workflow long_read_rna_quantification_kallisto.nf [required] [optional] [--help]

    required arguments:
        -c                                      :   Nextflow .config file.
        -w                                      :   Nextflow work directory path.
        --samples_tsv_file                      :   TSV file with the following columns:
                                                    'sample_id', 'fastq_file'.
        --kallisto_index_file                   :   Kallisto index file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/kallisto/v0.51.1/gencode_v41/kallisto_gencode_v41_index_k63.idx).
        --gtf_file                              :   GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf).
        --t2g_file                              :   T2G file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/kallisto/v0.51.1/gencode_v41/gencode_v41.t2g).
        --params_kallisto_bus                   :   Kallisto bus parameters (default: '"-x bulk --threshold 0.8"').
        --params_bustools_sort                  :   Bustools sort parameters (default: '""').
        --params_bustools_count                 :   Bustools count parameters (default: '"--cm -m"').
                                                    Note that the parameters need to be wrapped in quotes
                                                    and a space at the end of the string is necessary.
        --params_kallisto_quanttcc              :   Kallisto quant-tcc parameters (default: '"-P PacBio --matrix-to-files"').
        --output_dir                            :   Directory to which output files will be copied.

    optional arguments:
        --delete_work_dir                       :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                        :   ${params.samples_tsv_file}
        kallisto_index_file                     :   ${params.kallisto_index_file}
        gtf_file                                :   ${params.gtf_file}
        t2g_file                                :   ${params.t2g_file}
        params_kallisto_bus                     :   ${params_kallisto_bus}
        params_bustools_sort                    :   ${params_bustools_sort}
        params_bustools_count                   :   ${params_bustools_count}
        params_kallisto_quanttcc                :   ${params_kallisto_quanttcc}
        output_dir                              :   ${params.output_dir}
        delete_work_dir                         :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.fastq_file}") }
    .set { input_fastq_files_ch }

// Step 5. Workflow
workflow LONG_READ_RNA_QUANTIFICATION_KALLISTO {
    take:
        input_fastq_files_ch
        kallisto_index_file
        gtf_file
        t2g_file
        params_kallisto_bus
        params_bustools_sort
        params_bustools_count
        params_kallisto_quanttcc
        output_dir
    main:
        runKallistoQuantTccLongReads(
            input_fastq_files_ch,
            kallisto_index_file,
            gtf_file,
            t2g_file,
            params_kallisto_bus,
            params_bustools_sort,
            params_bustools_count,
            params_kallisto_quanttcc,
            output_dir
        )
    emit:
        runKallistoQuantTccLongReads.out.f
}

workflow {
    LONG_READ_RNA_QUANTIFICATION_KALLISTO(
        input_fastq_files_ch,
        params.kallisto_index_file,
        params.gtf_file,
        params.t2g_file,
        params_kallisto_bus,
        params_bustools_sort,
        params_bustools_count,
        params_kallisto_quanttcc,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
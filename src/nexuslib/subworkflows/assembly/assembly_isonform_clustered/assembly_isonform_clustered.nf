#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runIsonFormClustered }    from '../../../tools/isonform'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
// Reproduces isON_pipeline.sh's 'pacbio' mode, which invokes isONform_parallel
// with exactly these settings (see /isON_pipeline.sh in the container). The
// clustering half of that pipeline is replaced here by the supplied cluster TSV.
params.params_isonform      = '--k 20 --w 31 --xmin 14 --xmax 80 --exact_instance_limit 50 --max_seqs_to_spoa 200 --delta_len 10 --delta_iso_len_3 30 --delta_iso_len_5 50 --iso_abundance 3'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ASSEMBLY_ISONFORM_CLUSTERED {
    take:
        input_files_ch            // channel: [val(sample_id), path(fastq_file), path(tsv_file), val(cluster_method)]
        params_isonform
        output_dir

    main:
        runIsonFormClustered(
            input_files_ch,
            params_isonform,
            output_dir
        )

    emit:
        runIsonFormClustered.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===========================================================
             Assemble clustered long-read RNA FASTQ files using isONform
             ===========================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Assemble clustered transcripts using a custom batch script for isONform.

        usage: nexus run --nf-workflow assembly_isonform_clustered.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file', 'tsv_file' (columns: 'cluster_id', 'read_name'),
                                        'cluster_method' (name of the tool that produced tsv_file, e.g.
                                        'isonclust3'; keeps repeat clusterings of one sample apart).
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_isonform       :   isONform_parallel parameters (default: '"--k 20 --w 31 --xmin 14 --xmax 80 --exact_instance_limit 50 --max_seqs_to_spoa 200 --delta_len 10 --delta_iso_len_3 30 --delta_iso_len_5 50 --iso_abundance 3"').
                                        Note that the parameters need to be wrapped in quotes.
                                        Do NOT pass --write_fastq; the batch script reads back transcriptome.fasta.
                                        --t, --fastq_folder, --outfolder, --tmpdir and --split_wrt_batches are set
                                        by the batch script and must not be overridden here.

        """.stripIndent()
        exit 0
    }

    def params_isonform     = (params.params_isonform == true) ? '' : params.params_isonform

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_isonform         :   ${params_isonform}
    """.stripIndent()

    def cluster_rows = file(params.samples_tsv_file).splitCsv(header: true, sep: '\t')
    def missing_cm = cluster_rows.findAll { !(it.cluster_method ?: '').trim() }
    if (missing_cm) {
        error "ERROR: samples_tsv_file needs a non-blank 'cluster_method' column. " +
              "Missing for sample_id(s): ${missing_cm.collect { it.sample_id }.unique().join(', ')}"
    }
    def dup_cm = cluster_rows.groupBy { [it.sample_id, it.cluster_method] }
                             .findAll { k, v -> v.size() > 1 }
                             .keySet()
    if (dup_cm) {
        error "ERROR: duplicate (sample_id, cluster_method) rows would overwrite each " +
              "other: ${dup_cm.collect { it.join('/') }.join(', ')}"
    }

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}",
            "${row.tsv_file}",
            "${row.cluster_method}") }
        .set { input_files_ch }

    ASSEMBLY_ISONFORM_CLUSTERED(
        input_files_ch,
        params_isonform,
        params.output_dir
    )
}

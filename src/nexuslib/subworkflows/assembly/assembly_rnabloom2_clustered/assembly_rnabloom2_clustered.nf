#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runRnaBloom2LongReadClustered }    from '../../../tools/rnabloom2'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''

// Optional arguments
params.params_rnabloom2     = '--qual 20 --qual-avg 20 --mincov 3 -ntcard -savebf -chimera -lrpb'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ASSEMBLY_RNABLOOM2_CLUSTERED {
    take:
        input_files_ch            // channel: [val(sample_id), path(fastq_file), path(tsv_file), val(cluster_method)]
        params_rnabloom2
        output_dir

    main:
        runRnaBloom2LongReadClustered(
            input_files_ch,
            params_rnabloom2,
            output_dir
        )

    emit:
        runRnaBloom2LongReadClustered.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ============================================================
             Assemble clusterd long-read RNA FASTQ files using RNA-Bloom2
             ============================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Assemble clustered transcripts using a custom batch script for RNA-Bloom2.

        usage: nexus run --nf-workflow assembly_rnabloom2_clustered.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fastq_file', 'tsv_file' (columns: 'cluster_id', 'read_name'),
                                        'cluster_method' (name of the tool that produced tsv_file, e.g.
                                        'isonclust3'; keeps repeat clusterings of one sample apart).
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_rnabloom2      :   RNA-Bloom2 parameters (default: '"--qual 20 --qual-avg 20 --mincov 3 -ntcard -savebf -chimera -lrpb"').
                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_rnabloom2    = (params.params_rnabloom2 == true) ? '' : params.params_rnabloom2

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        params_rnabloom2        :   ${params_rnabloom2}
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

    ASSEMBLY_RNABLOOM2_CLUSTERED(
        input_files_ch,
        params_rnabloom2,
        params.output_dir
    )
}

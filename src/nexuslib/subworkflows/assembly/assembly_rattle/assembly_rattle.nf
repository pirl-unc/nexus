#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runRattleCluster }    from '../../../tools/rattle'
include { runRattleCorrect }    from '../../../tools/rattle'
include { runRattlePolish }     from '../../../tools/rattle'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                     = ''

// Required arguments
params.samples_tsv_file         = ''
params.output_dir               = ''

// Optional arguments
params.params_rattle_cluster    = '--iso --rna'
params.params_rattle_correct    = ''
params.params_rattle_polish     = '--rna'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ASSEMBLY_RATTLE {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        params_rattle_cluster
        params_rattle_correct
        params_rattle_polish
        output_dir

    main:
        // Step 1. Cluster
        runRattleCluster(
            input_fastq_files_ch,
            params_rattle_cluster,
            output_dir
        )

        // Step 2. Correct
        // Join original FASTQ with cluster output by sample_id
        // to produce: [sample_id, fastq_file, clusters.out]
        correct_input_ch = input_fastq_files_ch
            .join(runRattleCluster.out.f)
        runRattleCorrect(
            correct_input_ch,
            params_rattle_correct,
            output_dir
        )

        // Step 3. Polish
        runRattlePolish(
            runRattleCorrect.out.consensi,
            params_rattle_polish,
            output_dir
        )

    emit:
        runRattlePolish.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===============================================
             Assemble long-read RNA FASTQ files using RATTLE
             ===============================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Cluster reads using RATTLE.
            2. Correct reads using RATTLE.
            3. Polish consensus sequences using RATTLE.

        usage: nexus run --nf-workflow assembly_rattle.nf [required] [optional] [--help]

        required arguments:
            -c                              :   Nextflow .config file.
            -w                              :   Nextflow work directory path.
            --samples_tsv_file              :   TSV file with the following columns:
                                                'sample_id', 'fastq_file'.
            --output_dir                    :   Directory to which output files will be copied.

        optional arguments:
            --params_rattle_cluster         :   RATTLE cluster parameters (default: '"--iso --rna"').
                                                Note that the parameters need to be wrapped in quotes.
            --params_rattle_correct         :   RATTLE correct parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
            --params_rattle_polish          :   RATTLE polish parameters (default: '"--rna"').
                                                Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_rattle_cluster = (params.params_rattle_cluster == true) ? '' : params.params_rattle_cluster
    def params_rattle_correct = (params.params_rattle_correct == true) ? '' : params.params_rattle_correct
    def params_rattle_polish  = (params.params_rattle_polish == true) ? '' : params.params_rattle_polish

    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        params_rattle_cluster           :   ${params_rattle_cluster}
        params_rattle_correct           :   ${params_rattle_correct}
        params_rattle_polish            :   ${params_rattle_polish}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_fastq_files_ch }

    ASSEMBLY_RATTLE(
        input_fastq_files_ch,
        params_rattle_cluster,
        params_rattle_correct,
        params_rattle_polish,
        params.output_dir
    )
}

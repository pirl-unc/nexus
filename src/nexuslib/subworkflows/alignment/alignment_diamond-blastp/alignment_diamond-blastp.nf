#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runDiamondBlastp }    from '../../../tools/diamond'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                     = ''

// Required arguments
params.samples_tsv_file         = ''
params.output_dir               = ''
params.database_dmnd_file       = ''

// Optional arguments
params.params_diamond_blastp    = ''

if (params.params_diamond_blastp == true) {
    params_diamond_blastp = ''
} else {
    params_diamond_blastp = params.params_diamond_blastp
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ==============================================
         Align peptide fasta files using Diamond Blastp
         ==============================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Blast peptide sequences (fasta or fasta.gz files) using Diamond Blastp.

    usage: nexus run --nf-workflow alignment_diamond-blastp.nf [required] [optional] [--help]

    required arguments:
        -c                          :   Nextflow .config file.
        -w                          :   Nextflow work directory path.
        --samples_tsv_file          :   TSV file with the following columns:
                                        'sample_id', 'fasta_file'.
        --output_dir                :   Directory to which output files will be copied.
        --database_dmnd_file        :   Diamond-prepared local Blast database file.

    optional arguments:
        --params_diamond_blastp     :   Diamond blastp parameters (default: '""').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file            :   ${params.samples_tsv_file}
        output_dir                  :   ${params.output_dir}
        database_dmnd_file          :   ${params.database_dmnd_file}
        params_diamond_blastp       :   ${params_diamond_blastp}
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
workflow ALIGNMENT_DIAMOND_BLASTP {
    take:
        input_fasta_files_ch            // channel: [val(sample_id), path(fasta_file)]
        database_dmnd_file
        params_diamond_blastp
        output_dir

    main:
        runDiamondBlastp(
            input_fasta_files_ch,
            database_dmnd_file,
            params_diamond_blastp,
            output_dir
        )

    emit:
        runDiamondBlastp.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    ALIGNMENT_DIAMOND_BLASTP(
        input_fasta_files_ch,
        params.database_dmnd_file,
        params_diamond_blastp,
        params.output_dir
    )
}

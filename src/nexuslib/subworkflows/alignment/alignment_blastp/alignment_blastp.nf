#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runBlastp }   from '../../../tools/blastp'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.output_dir           = ''
params.blastdb_dir          = ''

// Optional arguments
params.blastdb_name         = 'nr'
params.params_blastp        = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ALIGNMENT_BLASTP {
    take:
        input_fasta_files_ch            // channel: [val(sample_id), path(fasta_file)]
        blastdb_dir
        blastdb_name
        params_blastp
        output_dir

    main:
        runBlastp(
            input_fasta_files_ch,
            blastdb_dir,
            blastdb_name,
            params_blastp,
            output_dir
        )

    emit:
        runBlastp.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ======================================
             Align peptide fasta files using Blastp
             ======================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Blast peptide sequences (fasta or fasta.gz files) using Blastp.

        usage: nexus run --nf-workflow alignment_blastp.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id', 'fasta_file'.
            --output_dir            :   Directory to which output files will be copied.
            --blastdb_dir           :   Local Blast database path. The database can be downloaded by the
                                        following command: update_blastdb.pl --source gcp --decompress nr --num_threads 16

        optional arguments:
            --blastdb_name          :   Local Blast database name (default: 'nr').
            --params_blastp         :   Blastp parameters (default: '""').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
        """.stripIndent()
        exit 0
    }

    def params_blastp = (params.params_blastp == true) ? '' : params.params_blastp

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        blastdb_dir             :   ${params.blastdb_dir}
        blastdb_name            :   ${params.blastdb_name}
        params_blastp           :   ${params_blastp}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fasta_file}") }
        .set { input_fasta_files_ch }

    ALIGNMENT_BLASTP(
        input_fasta_files_ch,
        params.blastdb_dir,
        params.blastdb_name,
        params_blastp,
        params.output_dir
    )
}

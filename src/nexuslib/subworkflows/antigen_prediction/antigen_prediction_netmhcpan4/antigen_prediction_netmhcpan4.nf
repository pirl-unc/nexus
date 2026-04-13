#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runNetMHCpan4 }   from '../../../tools/netmhcpan4'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.netmhcpan4_home_dir   = ''
params.netmhcpan4            = ''
params.samples_tsv_file      = ''
params.output_dir            = ''

// Optional arguments
params.params_netmhcpan4     = '-BA -s -l 8,9,10,11'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ANTIGEN_PREDICTION_NETMHCPAN4 {
    take:
        input_netmhcpan_files_ch             // channel: [val(sample_id), path(fasta_file), path(mhc_alleles)]
        netmhcpan4_home_dir
        netmhcpan4
        params_netmhcpan4
        output_dir

    main:
        runNetMHCpan4(
            input_netmhcpan_files_ch,
            netmhcpan4_home_dir,
            netmhcpan4,
            params_netmhcpan4,
            output_dir
        )

    emit:
        runNetMHCpan4.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ================================================
             Predict antigen presentation using NetMHCpan 4.x
             ================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run netMHCpan command.

        usage: nexus run --nf-workflow antigen_prediction_netmhcpan4.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns:
                                        'sample_id',
                                        'fasta_file',
                                        'mhc_alleles'.
            --output_dir            :   Directory to which output files will be copied.
            --netmhcpan4_home_dir   :   NetMHCpan home directory path.
            --netmhcpan4            :   NetMHCpan executable path.

        optional arguments:
            --params_netmhcpan4     :   netmhcpan parameters (default: '""').
                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_netmhcpan4 = (params.params_netmhcpan4 == true) ? '' : params.params_netmhcpan4

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        output_dir              :   ${params.output_dir}
        netmhcpan4_home_dir     :   ${params.netmhcpan4_home_dir}
        netmhcpan4              :   ${params.netmhcpan4}
        params_netmhcpan4       :   ${params_netmhcpan4}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fasta_file}",
            "${row.mhc_alleles}") }
        .set { input_netmhcpan_files_ch }

    ANTIGEN_PREDICTION_NETMHCPAN4(
        input_netmhcpan_files_ch,
        params.netmhcpan4_home_dir,
        params.netmhcpan4,
        params_netmhcpan4,
        params.output_dir
    )
}

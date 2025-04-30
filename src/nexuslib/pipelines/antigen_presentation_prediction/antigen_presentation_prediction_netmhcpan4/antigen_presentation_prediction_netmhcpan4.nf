#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runNetMHCpan } from '../../modules/netmhcpan'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.netmhcpan_home_dir = ''
params.netmhcpan = ''
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.params_netmhcpan = '-BA -s -l 8,9,10,11'
params.delete_work_dir = false

if (params.params_netmhcpan == true) {
    params_netmhcpan = ''
} else {
    params_netmhcpan = params.params_netmhcpan
}

// Step 3. Print inputs and help
log.info """\
         ================================================
         Predict antigen presentation using NetMHCpan 4.x
         ================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run netMHCpan command.

    usage: nexus run --nf-workflow antigen_presentation_prediction_netmhcpan4.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --netmhcpan_home_dir                :   NetMHCpan home directory path.
        --netmhcpan                         :   NetMHCpan executable path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id',
                                                'fasta_file',
                                                'mhc_alleles'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --params_netmhcpan                  :   netmhcpan parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        netmhcpan_home_dir                  :   ${params.netmhcpan_home_dir}
        netmhcpan                           :   ${params.netmhcpan}
        params_netmhcpan                    :   ${params_netmhcpan}
        delete_work_dir                     :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.fasta_file}",
        "${row.mhc_alleles}") }
    .set { input_netmhcpan_files_ch }

// Step 5. Workflow
workflow ANTIGEN_PRESENTATION_PREDICTION_NETMHCPAN4 {
    take:
        input_netmhcpan_files_ch             // channel: [val(sample_id), path(fasta_file), path(mhc_alleles)]
        netmhcpan_home_dir
        netmhcpan
        params_netmhcpan
        output_dir

    main:
        runNetMHCpan(
            input_netmhcpan_files_ch,
            netmhcpan_home_dir,
            netmhcpan,
            params_netmhcpan,
            output_dir
        )
}

workflow {
    ANTIGEN_PRESENTATION_PREDICTION_NETMHCPAN4(
        input_netmhcpan_files_ch,
        params.netmhcpan_home_dir,
        params.netmhcpan,
        params_netmhcpan,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }       from '../../../tools/samtools'
include { runMantaGermline }            from '../../../tools/manta'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''

// Optional arguments
params.params_manta_config              = ''
params.params_manta_run                 = ''

if (params.params_manta_config == true) {
    params_manta_config = ''
} else {
    params_manta_config = params.params_manta_config
}
if (params.params_manta_run == true) {
    params_manta_run = ''
} else {
    params_manta_run = params.params_manta_run
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ==================================================================================
         Identify germline variants in paired-end read DNA sequencing BAM files using Manta
         ==================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run manta (tumor and normal mode).

    usage: nexus run --nf-workflow variant_calling_manta-germline.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id',
                                                'bam_file',
                                                'bam_bai_file'
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.

    optional arguments:
        --params_manta_config               :   Manta configManta.py parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
        --params_manta_run                  :   Manta runWorkflow.py parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_manta_config                 :   ${params_manta_config}
        params_manta_run                    :   ${params_manta_run}
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
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_MANTA_GERMLINE {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        params_manta_config
        params_manta_run
        output_dir

    main:
        runSamtoolsFaidxFasta(reference_genome_fasta_file)
        fasta_file      = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file  = runSamtoolsFaidxFasta.out.fai_file

        runMantaGermline(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            params_manta_config,
            params_manta_run,
            output_dir
        )

    emit:
        runMantaGermline.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_MANTA_GERMLINE(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params_manta_config,
        params_manta_run,
        params.output_dir
    )
}

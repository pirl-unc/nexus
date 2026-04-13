#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }   from '../../../tools/samtools'
include { runNanoCaller }           from '../../../tools/nanocaller'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir                   = ''
params.reference_genome_fasta_file  = ''

// Optional arguments
params.params_nanocaller            = '--preset ont --mode all --sequencing ont --phase --enable_whatshap'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_NANOCALLER {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        params_nanocaller
        output_dir

    main:
        runSamtoolsFaidxFasta(reference_genome_fasta_file)
        fasta_file                  = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file              = runSamtoolsFaidxFasta.out.fai_file

        runNanoCaller(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            params_nanocaller,
            output_dir
        )

    emit:
        runNanoCaller.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===========================================================================================
             Identify germline small DNA variants in long-read DNA sequencing BAM files using NanoCaller
             ===========================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run Nanomonsv.

        usage: nexus run --nf-workflow variant_calling_nanocaller.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id',
                                                    'bam_file',
                                                    'bam_bai_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.

        optional arguments:
            --params_nanocaller                 :   NanoCaller parameters (default: '"--preset ont --mode all --sequencing ont --phase --enable_whatshap"').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_nanocaller = (params.params_nanocaller == true) ? '' : params.params_nanocaller

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_nanocaller                   :   ${params_nanocaller}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    VARIANT_CALLING_NANOCALLER(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params_nanocaller,
        params.output_dir
    )
}

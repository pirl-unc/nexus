#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }        from '../../../tools/samtools'
include { runSvimAlignmentMode }    from '../../../tools/svim'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''

// Optional arguments
params.params_svim                      = '--min_mapq 20 --min_sv_size 30 --insertion_sequences --read_names --zmws'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_SVIM {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        params_svim
        output_dir

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file          = runSamtoolsFaidx.out.fasta
        fasta_fai_file      = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file      = runSamtoolsFaidx.out.gzi_file

        runSvimAlignmentMode(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            params_svim,
            output_dir
        )

    emit:
        runSvimAlignmentMode.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =============================================================================
             Identify structural variants in long-read DNA sequencing BAM files using SVIM
             =============================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run SVIM.

        usage: nexus run --nf-workflow variant_calling_svim.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'bam_file', 'bam_bai_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.

        optional arguments:
            --params_svim                       :   SVIM parameters (default: '"--min_mapq 20 --min_sv_size 30 --insertion_sequences --read_names --zmws"').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_svim = (params.params_svim == true) ? '' : params.params_svim

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_svim                         :   ${params_svim}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.bam_file}",
            "${row.bam_bai_file}") }
        .set { input_bam_files_ch }

    VARIANT_CALLING_SVIM(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params_svim,
        params.output_dir
    )
}

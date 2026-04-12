#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }   from '../../../tools/samtools'
include { runColorSV }              from '../../../tools/colorsv'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir                   = ''
params.reference_genome_fasta_file  = ''
params.filter_bed_file              = ''

// Optional arguments
params.params_colorsv_preprocess    = '--read-sep /'
params.params_colorsv_call          = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_COLORSV {
    take:
        input_gfa_files_ch             // channel: [val(sample_id), path(gfa_file), val(tumor_ids)]
        reference_genome_fasta_file
        filter_bed_file
        params_colorsv_preprocess
        params_colorsv_call
        output_dir

    main:
        runSamtoolsFaidxFasta(reference_genome_fasta_file)
        fasta_file      = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file  = runSamtoolsFaidxFasta.out.fai_file

        runColorSV(
            input_gfa_files_ch,
            fasta_file,
            fasta_fai_file,
            filter_bed_file,
            params_colorsv_preprocess,
            params_colorsv_call,
            output_dir
        )

    emit:
        runColorSV.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ========================================================================================
             Identify somatic structural variants in long-read DNA sequencing GFA files using colorSV
             ========================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run colorSV.

        usage: nexus run --nf-workflow variant_calling_colorsv.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'gfa_file', 'tumor_ids'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --filter_bed_file                   :   BED file of regions to exclude.

        optional arguments:
            --params_colorsv_preprocess         :   colorSV preprocess parameters (default: '"--read-sep /"').
                                                    Note that the parameters need to be wrapped in quotes.
            --params_colorsv_call               :   colorSV call parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_colorsv_preprocess = (params.params_colorsv_preprocess == true) ? '' : params.params_colorsv_preprocess
    def params_colorsv_call       = (params.params_colorsv_call == true) ? '' : params.params_colorsv_call

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        filter_bed_file                     :   ${params.filter_bed_file}
        params_colorsv_preprocess           :   ${params_colorsv_preprocess}
        params_colorsv_call                 :   ${params_colorsv_call}
    """.stripIndent()

    Channel
        .fromPath(params.samples_tsv_file)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.sample_id,
                row.gfa_file,
                row.tumor_ids
            )
        }
        .set { input_gfa_files_ch }

    input_gfa_files_ch.view()

    VARIANT_CALLING_COLORSV(
        input_gfa_files_ch,
        params.reference_genome_fasta_file,
        params.filter_bed_file,
        params_colorsv_preprocess,
        params_colorsv_call,
        params.output_dir
    )
}

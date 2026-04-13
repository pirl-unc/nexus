#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }                  from '../../../tools/samtools'
include { runOctopusSomatic }                      from '../../../tools/octopus'
include { decompressFile as decompressFasta }      from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.regions_txt_file                 = ''

// Optional arguments
params.params_octopus                   = '--min-mapping-quality 20 --min-supporting-reads 3'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_OCTOPUS_SOMATIC {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)]
        reference_genome_fasta_file
        regions_txt_file
        params_octopus
        output_dir

    main:
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidxFasta(decompressFasta.out.f)
        fasta_file                  = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file              = runSamtoolsFaidxFasta.out.fai_file

        runOctopusSomatic(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            regions_txt_file,
            params_octopus,
            output_dir
        )

    emit:
        runOctopusSomatic.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===================================================================================
             Identify somatic variants in paired-end read DNA sequencing BAM files using Octopus
             ===================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1.  Run Octopus.

        usage: nexus run --nf-workflow variant_calling_octopus-somatic.nf [required] [optional] [--help]

        required arguments:
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'.
            --output_dir                        :   Directory to which output files will be symlinked.
            --reference_genome_fasta_file       :   Reference genome FASTA file.

        optional arguments:
            --params_octopus                    :   Octopus somatic mode parameters (default: '"--min-mapping-quality 20 --min-supporting-reads 3"').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_octopus = (params.params_octopus == true) ? '' : params.params_octopus

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_octopus                      :   ${params_octopus}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.tumor_bam_file}",
            "${row.tumor_bam_bai_file}",
            "${row.normal_bam_file}",
            "${row.normal_bam_bai_file}") }
        .set { input_bam_files_ch }

    VARIANT_CALLING_OCTOPUS_SOMATIC(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.regions_txt_file,
        params_octopus,
        params.output_dir
    )
}

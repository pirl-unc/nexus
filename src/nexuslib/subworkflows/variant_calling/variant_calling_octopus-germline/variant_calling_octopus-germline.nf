#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidxFasta }                  from '../../../tools/samtools'
include { runOctopusGermline }                     from '../../../tools/octopus'
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

if (params.params_octopus == true) {
    params_octopus = ''
} else {
    params_octopus = params.params_octopus
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ====================================================================================
         Identify germline variants in paired-end read DNA sequencing BAM files using Octopus
         ====================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1.  Run Octopus.

    usage: nexus run --nf-workflow variant_calling_octopus-germline.nf [required] [optional] [--help]

    required arguments:
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be symlinked.
        --reference_genome_fasta_file       :   Reference genome FASTA file.

    optional arguments:
        --params_octopus                    :   Octopus germline mode parameters (default: '"--min-mapping-quality 20 --min-supporting-reads 3"').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_octopus                      :   ${params_octopus}
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
workflow VARIANT_CALLING_OCTOPUS_GERMLINE {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        regions_txt_file
        params_octopus
        output_dir

    main:
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidxFasta(decompressFasta.out.f)
        fasta_file                  = runSamtoolsFaidxFasta.out.fasta
        fasta_fai_file              = runSamtoolsFaidxFasta.out.fai_file

        runOctopusGermline(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            regions_txt_file,
            params_octopus,
            output_dir
        )

    emit:
        runOctopusGermline.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_OCTOPUS_GERMLINE(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.regions_txt_file,
        params_octopus,
        params.output_dir
    )
}

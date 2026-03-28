#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }    from '../../../tools/samtools'
include { runSniffles2 }        from '../../../tools/sniffles2'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''

// Optional arguments
params.params_sniffles2                 = '--minsupport 3 --minsvlen 30 --mapq 20 --output-rnames'

if (params.params_sniffles2 == true) {
    params_sniffles2 = ''
} else {
    params_sniffles2 = params.params_sniffles2
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ==================================================================================
         Identify structural variants in long-read DNA sequencing BAM files using Sniffles2
         ==================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run Sniffles2.

    usage: nexus run --nf-workflow variant_calling_sniffles2.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.

    optional arguments:
        --params_sniffles2                  :   Sniffles2 parameters (default: '"--minsupport 3 --minsvlen 30 --mapq 20 --output-rnames"').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        params_sniffles2                    :   ${params_sniffles2}
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
workflow VARIANT_CALLING_SNIFFLES2 {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        params_sniffles2
        output_dir

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file                  = runSamtoolsFaidx.out.fasta
        fasta_fai_file              = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file              = runSamtoolsFaidx.out.gzi_file

        runSniffles2(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            params_sniffles2,
            output_dir
        )

    emit:
        runSniffles2.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_SNIFFLES2(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params_sniffles2,
        params.output_dir
    )
}

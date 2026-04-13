#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }                        from '../../../tools/samtools'
include { prepareGencodeGtfFileForReditools2 }      from '../../../tools/reditools2'
include { runReditools2 }                           from '../../../tools/reditools2'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.reference_genes_gtf_file         = ''

// Optional arguments
params.params_reditools2                = ''
params.params_reditools_annotatetable   = '-s 4 -c 1,2,3 -n gencode'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_REDITOOLS2 {
    take:
        input_bam_files_ch              // channel: [val(sample_id), path(rna_bam_file), path(rna_bam_bai_file), path(dna_bam_file), path(dna_bam_bai_file)]
        reference_genome_fasta_file
        reference_genes_gtf_file
        params_reditools2
        params_reditools_annotatetable
        output_dir

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file          = runSamtoolsFaidx.out.fasta
        fasta_fai_file      = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file      = runSamtoolsFaidx.out.gzi_file

        prepareGencodeGtfFileForReditools2(reference_genes_gtf_file)
        gtf_file            = prepareGencodeGtfFileForReditools2.out.gtf_file
        tbi_file            = prepareGencodeGtfFileForReditools2.out.tbi_file

        runReditools2(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            fasta_gzi_file,
            gtf_file,
            params_reditools2,
            params_reditools_annotatetable,
            output_dir
        )

    emit:
        runReditools2.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ====================================================================================
             Identify A-to-I RNA editing in paired-read RNA sequencing BAM files using Reditools2
             ====================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run REDItools2 reditools.py for RNA.
            2. Run REDItools2 reditools.py for DNA.
            3. Run annotate_with_DNA.py (identify RNA-specific A-to-I RNA editing events).
            4. Annotate variants.

        usage: nexus run --nf-workflow variant_calling_reditools2.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id',
                                                    'rna_bam_file',
                                                    'rna_bam_bai_file',
                                                    'dna_bam_file',
                                                    'dna_bam_bai_file'
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --reference_genes_gtf_file          :   Reference genes GTF file.

        optional arguments:
            --params_reditools2                 :   reditools.py parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
            --params_reditools_annotatetable    :   Reditools AnnotateTable.py parameters (default: '"-s 4 -c 1,2,3 -n gencode"').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_reditools2              = (params.params_reditools2 == true) ? '' : params.params_reditools2
    def params_reditools_annotatetable = (params.params_reditools_annotatetable == true) ? '' : params.params_reditools_annotatetable

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_reditools2                   :   ${params_reditools2}
        params_reditools_annotatetable      :   ${params_reditools_annotatetable}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.rna_bam_file}",
            "${row.rna_bam_bai_file}",
            "${row.dna_bam_file}",
            "${row.dna_bam_bai_file}") }
        .set { input_bam_files_ch }

    VARIANT_CALLING_REDITOOLS2(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params_reditools2,
        params_reditools_annotatetable,
        params.output_dir
    )
}

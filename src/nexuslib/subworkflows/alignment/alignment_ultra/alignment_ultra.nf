#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }                       from '../../../tools/samtools'
include { runUltraIndex }                          from '../../../tools/ultra'
include { runUltra }                               from '../../../tools/ultra'
include { runSamtoolsCalmd }                       from '../../../tools/samtools'
include { runSamtoolsSort }                        from '../../../tools/samtools'
include { copyBamFile }                            from '../../../tools/utils'
include { decompressFile as decompressFasta }      from '../../../tools/utils'
include { decompressFile as decompressGtf }        from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                         = ''

// Required arguments
params.samples_tsv_file             = ''
params.output_dir                   = ''
params.reference_genome_fasta_file  = ''
params.reference_genes_gtf_file     = ''

// Optional arguments
params.params_ultra_index = '--disable_infer'
params.params_ultra = '--isoseq'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ALIGNMENT_ULTRA {
    take:
        input_fastq_files_ch            // channel: [val(sample_id), path(fastq_file)]
        reference_genome_fasta_file
        reference_genes_gtf_file
        params_ultra_index
        params_ultra
        output_dir

    main:
        // Step 1. Decompress reference files
        decompressFasta(reference_genome_fasta_file)
        decompressGtf(reference_genes_gtf_file)

        // Step 2. Index reference genome FASTA file
        runSamtoolsFaidx(decompressFasta.out.f)
        fasta_file          = runSamtoolsFaidx.out.fasta
        fasta_fai_file      = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file      = runSamtoolsFaidx.out.gzi_file

        // Step 3. Create uLTRA index
        runUltraIndex(
            decompressFasta.out.f,
            decompressGtf.out.f,
            params_ultra_index
        )

        // Step 3. Run uLTRA
        runUltra(
            input_fastq_files_ch,
            fasta_file,
            fasta_fai_file,
            runUltraIndex.out.f,
            params_ultra,
            output_dir
        )

        // Step 4. Run Samtools calmd
        runSamtoolsCalmd(
            runUltra.out.f,
            fasta_file,
            fasta_fai_file
        )

        // Step 5. Run Samtools sort
        runSamtoolsSort(runSamtoolsCalmd.out.f)

        // Step 6. Copy BAM files
        copyBamFile(
            runSamtoolsSort.out.f,
            output_dir
        )

    emit:
        runSamtoolsSort.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ======================================================
             Align long-read RNA sequencing fastq files using uLTRA
             ======================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Align reads (fastq.gz files) to a reference genome using uLTRA.
            2. Generate MD tags.
            3. Sort MD-tagged bam file.

        usage: nexus run --nf-workflow alignment_ultra.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id', 'fastq_file'.
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --reference_genes_gtf_file          :   Reference genes GTF file.

        optional arguments:
            --params_ultra_index                :   uLTRA index parameters (default: "--disable_infer").
            --params_ultra                      :   uLTRA align parameters (default: "--isoseq").
        """.stripIndent()
        exit 0
    }

    def params_ultra_index = (params.params_ultra_index == true) ? '' : params.params_ultra_index
    def params_ultra       = (params.params_ultra == true) ? '' : params.params_ultra

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file            :   ${params.reference_genes_gtf_file}
        params_ultra_index                  :   ${params_ultra_index}
        params_ultra                        :   ${params_ultra}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file}") }
        .set { input_fastq_files_ch }

    ALIGNMENT_ULTRA(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params_ultra_index,
        params_ultra,
        params.output_dir
    )
}

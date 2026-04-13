#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runStarIndex }                           from '../../../tools/star'
include { runStar }                                from '../../../tools/star'
include { decompressFile as decompressFasta }      from '../../../tools/utils'
include { decompressFile as decompressGtf }        from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help = ''

// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
params.reference_genome_fasta_file  = ''
params.reference_genes_gtf_file     = ''

// Optional arguments
params.params_star_genomegenerate   = '--genomeSAindexNbases 10 --sjdbOverhang 150'
params.params_star                  = '--genomeLoad NoSharedMemory --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic --chimSegmentMin 10 --chimOutType WithinBAM SoftClip Junctions --chimMultimapNmax 20 --chimOutJunctionFormat 0'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow ALIGNMENT_STAR {
    take:
        input_fastq_files_ch
        reference_genome_fasta_file
        reference_genes_gtf_file
        params_star_genomegenerate
        params_star
        output_dir

    main:
        decompressFasta(reference_genome_fasta_file)
        decompressGtf(reference_genes_gtf_file)

        runStarIndex(
            decompressFasta.out.f,
            decompressGtf.out.f,
            params_star_genomegenerate
        )
        runStar(
            input_fastq_files_ch,
            runStarIndex.out.f,
            params_star,
            output_dir
        )

    emit:
        runStar.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ======================================================
             Align paired-end RNA sequencing fastq files using STAR
             ======================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Index reference genome FASTA and GTF files.
            2. Align paired-end reads to a reference genome index using STAR.

        usage: nexus run --nf-workflow alignment_star.nf [required] [optional] [--help]

        required arguments:
            -c                              :   Nextflow .config file.
            -w                              :   Nextflow work directory path.
            --samples_tsv_file              :   TSV file with the following columns:
                                                'sample_id', 'fastq_file_1', 'fastq_file_2'.
            --output_dir                    :   Directory to which output files will be copied.
            --reference_genome_fasta_file   :   Reference genome FASTA file.
            --reference_genes_gtf_file      :   Reference genes GTF file.

        optional arguments:
            --params_star_genomegenerate    :   STAR genomeGenerate parameter (default: "--genomeSAindexNbases 10 --sjdbOverhang 150").
                                                Note that the parameters need to be wrapped in quotes.
            --params_star                   :   STAR parameters (default: "--genomeLoad NoSharedMemory --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic --chimSegmentMin 10 --chimOutType WithinBAM SoftClip Junctions --chimMultimapNmax 20 --chimOutJunctionFormat 0").
                                                Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_star_genomegenerate = (params.params_star_genomegenerate == true) ? '' : params.params_star_genomegenerate
    def params_star                = (params.params_star == true) ? '' : params.params_star

    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        reference_genome_fasta_file     :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file        :   ${params.reference_genes_gtf_file}
        params_star_genomegenerate      :   ${params_star_genomegenerate}
        params_star                     :   ${params_star}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.fastq_file_1}",
            "${row.fastq_file_2}") }
        .set { input_fastq_files_ch }

    ALIGNMENT_STAR(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params_star_genomegenerate,
        params_star,
        params.output_dir
    )
}

#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runStar } from '../../modules/star'
include { runArriba } from '../../modules/arriba'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.reference_genome_fasta_fai_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai'
params.star_index = '/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/star/hg38_149bp_overhang/'
params.gtf_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf'
params.protein_domains_gff3_file = '/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/arriba/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3'
params.params_arriba = '-S 3 -f blacklist -i chr*'
params.params_star = '--genomeLoad NoSharedMemory --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --twopassMode Basic --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50'
params.delete_work_dir = false

if (params.params_arriba == true) {
    params_arriba = ''
} else {
    params_arriba = params.params_arriba
}
if (params.params_star == true) {
    params_star = ''
} else {
    params_star = params.params_star
}

// Step 3. Print inputs and help
log.info """\
         ===============================================================================
         Identify fusion genes in paired-read RNA sequencing FASTQ.GZ files using Arriba
         ===============================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Align paired-end reads to a reference genome index using STAR.
        2. Run Arriba.

    usage: nexus run --nf-workflow paired-end_read_rna_variant_calling_arriba.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns:
                                                'sample_id', 'fastq_file_1', 'fastq_file_2'.
        --output_dir                        :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
        --star_index                        :   Reference genome STAR index (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/star/hg38_100bp_overhang/).
        --gtf_file                          :   Reference transcriptome GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf).
        --protein_domains_gff3_file         :   Protein domains GFF3 file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/arriba/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3).
        --params_arriba                     :   Arriba parameters (default: '"-S 3 -f blacklist -i chr*"').
                                                Note that the parameters need to be wrapped in quotes.
        --params_star                       :   STAR parameters (default: '"--genomeLoad NoSharedMemory --readFilesCommand zcat --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50"').
                                                Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                   :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        reference_genome_fasta_fai_file     :   ${params.reference_genome_fasta_fai_file}
        star_index                          :   ${params.star_index}
        gtf_file                            :   ${params.gtf_file}
        protein_domains_gff3_file           :   ${params.protein_domains_gff3_file}
        params_arriba                       :   ${params_arriba}
        params_star                         :   ${params_star}
        delete_work_dir                     :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.fastq_file_1}",
        "${row.fastq_file_2}") }
    .set { input_fastq_files_ch }

// Step 5. Workflow
workflow PAIRED_END_READ_RNA_VARIANT_CALLING_ARRIBA {
    take:
        input_fastq_files_ch             // channel: [val(sample_id), path(fastq_file_1), path(fastq_file_2)]
        reference_genome_fasta_file
        reference_genome_fasta_fai_file
        star_index
        gtf_file
        protein_domains_gff3_file
        params_arriba
        params_star
        output_dir

    main:
        runStar(
            input_fastq_files_ch,
            star_index,
            params_star,
            output_dir
        )
        runArriba(
            runStar.out.f,
            reference_genome_fasta_file,
            reference_genome_fasta_fai_file,
            gtf_file,
            protein_domains_gff3_file,
            params_arriba,
            output_dir
        )
}

workflow {
    PAIRED_END_READ_RNA_VARIANT_CALLING_ARRIBA(
        input_fastq_files_ch,
        params.reference_genome_fasta_file,
        params.reference_genome_fasta_fai_file,
        params.star_index,
        params.gtf_file,
        params.protein_domains_gff3_file,
        params_arriba,
        params_star,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}

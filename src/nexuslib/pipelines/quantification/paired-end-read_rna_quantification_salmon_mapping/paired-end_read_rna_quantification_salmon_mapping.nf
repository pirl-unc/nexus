#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runSalmonPairedEndMappingMode } from '../../modules/salmon'

// Step 2. Input parameters
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_transcripts_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.transcripts.fa'
params.gtf_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf.gz'
params.params_salmon_index = '--gencode'
params.params_salmon_quant = '--libType IU --seqBias --gcBias --posBias'
params.delete_work_dir = false

if (params.params_salmon_index == true) {
    params_salmon_index = ''
} else {
    params_salmon_index = params.params_salmon_index
}
if (params.params_salmon_quant == true) {
    params_salmon_quant = ''
} else {
    params_salmon_quant = params.params_salmon_quant
}

// Step 3. Print inputs and help
log.info """\
         ===================================================
         Quantify RNA in paired-end FASTQ files using Salmon
         ===================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Quantify RNA in paired-end FASTQ files using Salmon.

    usage: nexus run --nf-workflow paired-end_read_rna_quantification_salmon_mapping.nf [required] [optional] [--help]

    required arguments:
        -c                                      :   Nextflow .config file.
        -w                                      :   Nextflow work directory path.
        --samples_tsv_file                      :   TSV file with the following columns:
                                                    'sample_id', 'fastq_file_1', 'fastq_file_2'.
        --output_dir                            :   Directory to which output files will be copied.

    optional arguments:
        --reference_transcripts_fasta_file      :   Reference transcripts FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.transcripts.fa).
        --gtf_file                              :   GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf.gz).
        --params_salmon_index                   :   Salmon index parameters (default: '"--gencode "').
                                                    Note that the parameters need to be wrapped in quotes.
        --params_salmon_quant                   :   Salmon parameters (default: '"--libType IU --seqBias --gcBias --posBias"').
                                                    Note that the parameters need to be wrapped in quotes.
        --delete_work_dir                       :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                        :   ${params.samples_tsv_file}
        output_dir                              :   ${params.output_dir}
        reference_transcripts_fasta_file        :   ${params.reference_transcripts_fasta_file}
        gtf_file                                :   ${params.gtf_file}
        params_salmon_index                     :   ${params_salmon_index}
        params_salmon_quant                     :   ${params_salmon_quant}
        delete_work_dir                         :   ${params.delete_work_dir}
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
workflow PAIRED_END_RNA_QUANTIFICATION_SALMON_MAPPING {
    take:
        input_fastq_files_ch
        reference_transcripts_fasta_file
        gtf_file
        params_salmon_index
        params_salmon_quant
        output_dir
    main:
        runSalmonPairedEndMappingMode(
            input_fastq_files_ch,
            reference_transcripts_fasta_file,
            gtf_file,
            params_salmon_index,
            params_salmon_quant,
            output_dir
        )
    emit:
        runSalmonPairedEndMappingMode.out.f
}

workflow {
    PAIRED_END_RNA_QUANTIFICATION_SALMON_MAPPING(
        input_fastq_files_ch,
        params.reference_transcripts_fasta_file,
        params.gtf_file,
        params_salmon_index,
        params_salmon_quant,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
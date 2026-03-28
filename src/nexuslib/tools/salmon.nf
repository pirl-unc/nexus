#!/usr/bin/env nextflow

process runSalmonIndex {

    label 'salmon_index'
    debug true

    input:
        path(reference_transcripts_fasta_file)
        val(params_salmon_index)

    output:
        path("salmon_index/"), emit: f

    script:
        """
        mkdir -p salmon_index/
        salmon index \
            --transcripts $reference_transcripts_fasta_file \
            --index salmon_index/ \
            --threads ${task.cpus} \
            $params_salmon_index
        """
}

process runSalmonFastqMode {

    label 'salmon'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        path(reference_transcripts_fasta_file)
        path(salmon_index_dir)
        val(params_salmon_quant)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_salmon_output/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_salmon_output/
        salmon quant \
            $params_salmon_quant \
            --index $salmon_index_dir \
            --mates1 $fastq_file_1 \
            --mates2 $fastq_file_2 \
            --geneMap $reference_transcripts_fasta_file \
            --output ${sample_id}_salmon_output/ \
            --threads ${task.cpus}
        """
}

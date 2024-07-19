#!/usr/bin/env nextflow

process runSeq2HLA {

    label 'seq2hla'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(params_seq2hla)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}-Class*"), emit: f

    script:
        """
        seq2HLA \
          -1 $fastq_file_1 \
          -2 $fastq_file_2 \
          -r $sample_id \
          -p ${task.cpus} \
          $params_seq2hla
        """
}
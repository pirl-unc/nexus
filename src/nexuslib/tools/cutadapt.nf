#!/usr/bin/env nextflow

process runCutadapt {

    label 'cutadapt'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(params_cutadapt)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${sample_id}_R1.trim.fastq.gz"),
              path("${sample_id}_R2.trim.fastq.gz"),
              emit: f

    script:
        """
        cutadapt \
            $params_cutadapt \
            -j ${task.cpus} \
            -o ${sample_id}_R1.trim.fastq.gz \
            -p ${sample_id}_R2.trim.fastq.gz \
            ${fastq_file_1} ${fastq_file_2}
        """
}
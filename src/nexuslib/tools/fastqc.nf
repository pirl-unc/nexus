#!/usr/bin/env nextflow

process runFastQCSingleEndRead {

    label 'fastqc'
    tag "${sample_id}"
    debug true

   publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(output_dir)

    output:
        tuple val(sample_id), path("fastqc/"), emit: f

    script:
        """
        mkdir -p fastqc/
        fastqc \
            $fastq_file \
            --memory ${task.fastqc_memory.toMega()} \
            --threads ${task.cpus} \
            -o fastqc/
        """
}

process runFastQCPairedEndRead {

    label 'fastqc'
    tag "${sample_id}"
    debug true

   publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(output_dir)

    output:
        tuple val(sample_id), path("fastqc/"), emit: f

    script:
        """
        mkdir -p fastqc/
        fastqc \
            $fastq_file_1 $fastq_file_2 \
            --memory ${task.fastqc_memory.toMega()} \
            --threads ${task.cpus} \
            -o fastqc/
        """
}

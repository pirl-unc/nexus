#!/usr/bin/env nextflow

process runFastQCSingleEndRead {

    label 'fastqc'
    tag "${sample_id}"
    debug true

   publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(output_dir)

    output:
        tuple val(sample_id), path("*.html"), path("*.zip"), emit: f

    script:
        """
        mkdir -p output/
        fastqc \
            $fastq_file \
            --memory ${task.fastqc_memory.toMega()} \
            --threads ${task.cpus} \
            -o output/
        mv output/* .
        """
}

process runFastQCPairedEndRead {

    label 'fastqc'
    tag "${sample_id}"
    debug true

   publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(output_dir)

    output:
        tuple val(sample_id), path("*.html"), path("*.zip"), emit: f

    script:
        """
        mkdir -p output/
        fastqc \
            $fastq_file_1 $fastq_file_2 \
            --memory ${task.fastqc_memory.toMega()} \
            --threads ${task.cpus} \
            -o output/
        mv output/* .
        """
}

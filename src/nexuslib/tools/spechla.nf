#!/usr/bin/env nextflow

process runSpecHLAShortRead {

    label 'spechla'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(params_spechla)
        val(output_dir)

    output:
        tuple val(sample_id), path("spechla/"), emit: f

    script:
        """
        mkdir -p spechla/
        spechla \
            -n $sample_id \
            -1 $fastq_file_1 \
            -2 $fastq_file_2 \
            -o spechla/ \
            -j ${task.cpus} \
            $params_spechla
        """
}

process runSpecHLALongRead {

    label 'spechla'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_spechla)
        val(output_dir)

    output:
        tuple val(sample_id), path("spechla/"), emit: f

    script:
        """
        mkdir -p spechla/
        spechla-long-read \
            -n $sample_id \
            -r $fastq_file \
            -o spechla/ \
            -j ${task.cpus} \
            $params_spechla
        """
}

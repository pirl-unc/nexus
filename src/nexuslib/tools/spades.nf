#!/usr/bin/env nextflow

process runRNASpades {

    label 'spades'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(params_rnaspades)
        val(output_dir)

    output:
        tuple val(sample_id), path("rnaspades/"), emit: f

    script:
        """
        mkdir -p rnaspades/
        rnaspades.py \
            -1 $fastq_file_1 \
            -2 $fastq_file_2 \
            -o rnaspades/ \
            -t ${task.cpus} \
            -m ${task.memory.toGiga()} \
            $params_rnaspades
        """
}

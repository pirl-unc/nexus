#!/usr/bin/env nextflow

process runTrinity {

    label 'trinity'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(params_trinity)
        val(output_dir)

    output:
        tuple val(sample_id), path("trinity/"), emit: f

    script:
        """
        mkdir -p trinity/trinity/
        Trinity \
            --seqType fq \
            --left $fastq_file_1 \
            --right $fastq_file_2 \
            --output trinity/trinity/ \
            --CPU ${task.cpus} \
            --max_memory ${task.memory.toGiga()}G \
            $params_trinity
        """
}

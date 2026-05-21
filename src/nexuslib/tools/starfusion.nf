#!/usr/bin/env nextflow

process runStarFusion {

    label 'star_fusion'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        path(genome_lib_dir)
        val(params_starfusion)
        val(output_dir)

    output:
        tuple val(sample_id), path("starfusion/"), emit: f

    script:
        """
        mkdir -p starfusion/
        STAR-Fusion \
            --left_fq $fastq_file_1 \
            --right_fq $fastq_file_2 \
            --genome_lib_dir \${PWD}/${genome_lib_dir} \
            --CPU ${task.cpus} \
            --output_dir starfusion/ \
            $params_starfusion
        """
}

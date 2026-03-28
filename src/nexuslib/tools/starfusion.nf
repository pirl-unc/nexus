#!/usr/bin/env nextflow

process runStarFusion {

    label 'star_fusion'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        path(genome_lib_dir)
        val(params_starfusion)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_starfusion_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_starfusion_outputs/
        STAR-Fusion \
            --left_fq $fastq_file_1 \
            --right_fq $fastq_file_2 \
            --genome_lib_dir \${PWD}/${genome_lib_dir} \
            --CPU ${task.cpus} \
            --output_dir ${sample_id}_starfusion_outputs/ \
            $params_starfusion
        """
}

#!/usr/bin/env nextflow

process runGeluster {

    label 'geluster'
    tag "${sample_id}"
    debug true
    // Stage inputs as real copies, not symlinks.
    stageInMode 'copy'

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_geluster)
        val(output_dir)

    output:
        tuple val(sample_id), path("geluster/"), emit: f

    script:
        """
        gunzip -c $fastq_file > ${sample_id}_long_read.fastq
        GeLuster \
            --reads ${sample_id}_long_read.fastq \
            --threads ${task.cpus} \
            --output_dir \${PWD}/geluster \
            $params_geluster
        """
}

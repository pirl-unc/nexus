#!/usr/bin/env nextflow

process runStringTie2 {

    label 'stringtie2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file)
        path(gtf_file)
        val(params_stringtie2)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_stringtie2.gtf"), emit: f

    script:
        """
        stringtie \
            -o ${sample_id}_stringtie2.gtf \
            -G $gtf_file \
            -p ${task.cpus} \
            $params_stringtie2 \
            $bam_file
        """
}

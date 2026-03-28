#!/usr/bin/env nextflow

process runLonggf {

    label 'longgf'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file)
        path(gtf_file)
        val(params_longgf)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_longgf_results.txt"), emit: f

    script:
        """
        LongGF \
            $bam_file \
            $gtf_file \
            $params_longgf > ${sample_id}_longgf_results.txt
        """
}
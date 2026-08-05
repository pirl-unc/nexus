#!/usr/bin/env nextflow

process runFufihla {

    label 'fufihla'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_fufihla)
        val(output_dir)

    output:
        tuple val(sample_id), path("fufihla/"), emit: f

    script:
        """
        mkdir -p fufihla/
        fufihla \
            --fa $fastq_file \
            --out fufihla \
            $params_fufihla
        """
}

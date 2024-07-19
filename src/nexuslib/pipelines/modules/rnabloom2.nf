#!/usr/bin/env nextflow

process runRnaBloom2LongRead {

    label 'rnabloom2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_rnabloom2)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_rnabloom2_output/"), emit: f

    script:
        """
        mkdir ${sample_id}_rnabloom2_output/
        java -jar -Xmx${task.java_max_mem.toGiga()}G /opt/rnabloom2/RNA-Bloom.jar \
            -long $fastq_file \
            --threads ${task.cpus} \
            --outdir ${sample_id}_rnabloom2_output/ \
            $params_rnabloom2
        """
}
#!/usr/bin/env nextflow

process runRnaBloom2LongRead {

    label 'rnabloom2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_rnabloom2)
        val(output_dir)

    output:
        tuple val(sample_id), path("rnabloom2/"), emit: f

    script:
        """
        mkdir rnabloom2/
        java -jar -Xmx${task.java_max_mem.toGiga()}G /opt/rnabloom2/RNA-Bloom.jar \
            -long $fastq_file \
            --threads ${task.cpus} \
            --outdir rnabloom2/ \
            $params_rnabloom2
        """
}

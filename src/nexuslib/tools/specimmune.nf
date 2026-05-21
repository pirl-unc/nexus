#!/usr/bin/env nextflow

process runSpecImmune {

    label 'specimmune'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_specimmune)
        val(output_dir)

    output:
        tuple val(sample_id), path("specimmune/"), emit: f

    script:
        """
        mkdir -p specimmune/
        python /opt/SpecImmune-0.0.3/scripts/main.py \
            -r $fastq_file \
            -n $sample_id \
            -o specimmune/ \
            --db /opt/SpecImmune-0.0.3/db/ \
            -j ${task.cpus} \
            $params_specimmune
        """
}
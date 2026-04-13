#!/usr/bin/env nextflow

process runRattleCluster {

    label 'rattle'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_rattle_cluster)
        val(output_dir)

    output:
        tuple val(sample_id), path("clusters.out"), emit: f

    script:
        """
        rattle cluster \
            --input $fastq_file \
            --output . \
            -t ${task.cpus} \
            $params_rattle_cluster
        """
}

process runRattleCorrect {

    label 'rattle'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file), path(clusters_out_file)
        val(params_rattle_correct)
        val(output_dir)

    output:
        tuple val(sample_id), path("corrected.fq"), emit: corrected
        tuple val(sample_id), path("uncorrected.fq"), emit: uncorrected
        tuple val(sample_id), path("consensi.fq"), emit: consensi

    script:
        """
        rattle correct \
            --input $fastq_file \
            --clusters $clusters_out_file \
            --output . \
            -t ${task.cpus} \
            $params_rattle_correct
        """
}

process runRattlePolish {

    label 'rattle'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(consensi_file)
        val(params_rattle_polish)
        val(output_dir)

    output:
        tuple val(sample_id), path("transcriptome.fq.gz"), emit: f

    script:
        """
        rattle polish \
            --input $consensi_file \
            --output-folder . \
            -t ${task.cpus} \
            $params_rattle_polish
        gzip transcriptome.fq
        """
}

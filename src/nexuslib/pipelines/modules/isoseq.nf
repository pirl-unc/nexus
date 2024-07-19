#!/usr/bin/env nextflow

process runIsoseqCluster {

    label 'isoseq_cluster'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_isoseq-clustered.bam"), path("${sample_id}_isoseq-clustered.bam.bai"), emit: f

    script:
        """
        ls $bam_file > flnc.fofn
        isoseq cluster2 \
            --num-threads ${task.cpus} \
            flnc.fofn \
            ${sample_id}_isoseq-clustered.bam
        samtools index -b ${sample_id}_isoseq-clustered.bam ${sample_id}_isoseq-clustered.bam.bai
        """
}

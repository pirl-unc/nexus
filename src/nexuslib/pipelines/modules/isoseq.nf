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
        val(isoseq)
        val(output_dir)

    output:
        tuple val(sample_id), path({"${bam_file.simpleName}_isoseq-clustered.bam"}), emit: f

    script:
        """
        ls $bam_file > flnc.fofn
        $isoseq cluster2 \
            --num-threads ${task.cpus} \
            flnc.fofn \
            ${bam_file.simpleName}_isoseq-clustered.bam
        """
}

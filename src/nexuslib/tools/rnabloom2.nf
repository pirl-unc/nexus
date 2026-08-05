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

process runRnaBloom2LongReadClustered {

    label 'rnabloom2_clustered'
    tag "${sample_id}:${cluster_method}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file), path(tsv_file), val(cluster_method)
        val(params_rnabloom2)
        val(output_dir)

    output:
        tuple val(sample_id), val(cluster_method), path("${cluster_method}_rnabloom2_clustered/"), emit: f

    script:
        """
        mkdir -p ${cluster_method}_rnabloom2_clustered/
        python /opt/rnabloom2/run_rnabloom2_clustered.py \
            --jar /opt/rnabloom2/RNA-Bloom.jar \
            --fastq-file $fastq_file \
            --tsv-file $tsv_file \
            --output-dir \${PWD}/${cluster_method}_rnabloom2_clustered/ \
            --output-prefix ${sample_id}_${cluster_method} \
            --num-threads ${task.num_threads_per_worker} \
            --num-parallel ${task.num_parallel} \
            --xmx ${task.java_max_mem_per_worker.toGiga()}G \
            --extra-args="$params_rnabloom2"
        """
}

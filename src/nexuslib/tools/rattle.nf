#!/usr/bin/env nextflow

process runRattleCluster {

    label 'rattle'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/rattle/",
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
        path: "${output_dir}/${sample_id}/rattle/",
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
        path: "${output_dir}/${sample_id}/rattle/",
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
        if [ -s $consensi_file ]; then
            rattle polish \
                --input $consensi_file \
                --output-folder . \
                -t ${task.cpus} \
                $params_rattle_polish
        else
            echo "WARNING: 'rattle correct' produced no consensus sequences for ${sample_id} — emitting an empty transcriptome."
            : > transcriptome.fq
        fi
        gzip transcriptome.fq
        """
}

process runRattleClustered {

    label 'rattle_clustered'
    tag "${sample_id}:${cluster_method}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file), path(tsv_file), val(cluster_method)
        val(params_rattle_cluster)
        val(params_rattle_correct)
        val(params_rattle_polish)
        val(output_dir)

    output:
        tuple val(sample_id), val(cluster_method), path("${cluster_method}_rattle_clustered/"), emit: f

    script:
        """
        mkdir -p ${cluster_method}_rattle_clustered/
        python /run_rattle_clustered.py \
            --rattle rattle \
            --fastq-file $fastq_file \
            --tsv-file $tsv_file \
            --output-dir \${PWD}/${cluster_method}_rattle_clustered/ \
            --output-prefix ${sample_id}_${cluster_method} \
            --num-threads ${task.num_threads_per_worker} \
            --num-parallel ${task.num_parallel} \
            --cluster-extra-args="$params_rattle_cluster" \
            --correct-extra-args="$params_rattle_correct" \
            --polish-extra-args="$params_rattle_polish"
        """
}

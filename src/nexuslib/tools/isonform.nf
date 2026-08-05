#!/usr/bin/env nextflow

process runIsonForm {

    label 'isonform'
    tag "${sample_id}"
    debug true
    // Stage inputs as real copies, not symlinks.
    stageInMode 'copy'

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_ison_pipeline)
        val(output_dir)

    output:
        tuple val(sample_id), path("ison_pipeline/"), emit: f

    script:
        """
        mkdir -p ison_pipeline/
        gunzip -c $fastq_file > ${sample_id}_long_read.fastq
        bash /isON_pipeline.sh \
            --raw_reads \$PWD/${sample_id}_long_read.fastq \
            --outfolder \$PWD/ison_pipeline/ \
            --num_cores ${task.cpus} \
            $params_ison_pipeline
        """
}

process runIsonFormClustered {

    label 'isonform_clustered'
    tag "${sample_id}:${cluster_method}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file), path(tsv_file), val(cluster_method)
        val(params_isonform)
        val(output_dir)

    output:
        tuple val(sample_id), val(cluster_method), path("${cluster_method}_isonform_clustered/"), emit: f

    script:
        """
        mkdir -p ${cluster_method}_isonform_clustered/
        python /run_isonform_clustered.py \
            --isonform-parallel isONform_parallel \
            --fastq-file $fastq_file \
            --tsv-file $tsv_file \
            --output-dir \${PWD}/${cluster_method}_isonform_clustered/ \
            --output-prefix ${sample_id}_${cluster_method} \
            --num-parallel ${task.num_parallel} \
            --extra-args="$params_isonform"
        """
}

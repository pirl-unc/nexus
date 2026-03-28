#!/usr/bin/env nextflow

process runRmatsBamMode {

    label 'rmats'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genes_gtf_file)
        val(params_rmats)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_rmats_outputs/"), emit: f

    script:
        """
        mkdir -p temp/
        mkdir -p ${sample_id}_rmats_outputs/
        realpath $bam_file > path.txt
        rmats.py \
            --b1 path.txt \
            --gtf $reference_genes_gtf_file \
            --tmp \$PWD/temp/ \
            --nthread ${task.cpus} \
            --od ${sample_id}_rmats_outputs/ \
            $params_rmats
        """
}
#!/usr/bin/env nextflow

process runRmatsBamMode {

    label 'rmats'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genes_gtf_file)
        val(params_rmats)
        val(output_dir)

    output:
        tuple val(sample_id), path("rmats/"), emit: f

    script:
        """
        mkdir -p temp/
        mkdir -p rmats/
        realpath $bam_file > path.txt
        rmats.py \
            --b1 path.txt \
            --gtf $reference_genes_gtf_file \
            --tmp \$PWD/temp/ \
            --nthread ${task.cpus} \
            --od rmats/ \
            $params_rmats
        """
}
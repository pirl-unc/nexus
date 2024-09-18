#!/usr/bin/env nextflow

process runLiqaQuantify {

    label 'liqa'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(liqa_refgene_file)
        val(params_liqa_quantify)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_liqa_quantify_output"), emit: f

    script:
        """
        samtools view $bam_file -F 2308 -q 50 -O BAM -o ${bam_file.baseName}_filtered.bam
        liqa -task quantify \
            -refgene $liqa_refgene_file \
            -bam ${bam_file.baseName}_filtered.bam \
            -out ${sample_id}_liqa_quantify_output \
            $params_liqa_quantify
        """
}

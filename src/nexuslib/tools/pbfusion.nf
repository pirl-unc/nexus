#!/usr/bin/env nextflow

process runPbfusion {

    label 'pbfusion'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(gtf_file)
        val(params_pbfusion_discover)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_pbfusion*"), emit: f

    script:
        """
        mkdir -p outputs/
        pbfusion discover \
            -b $bam_file \
            -g $gtf_file \
            -o ${sample_id}_pbfusion \
            $params_pbfusion_discover
        """
}
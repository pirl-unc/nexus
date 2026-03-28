#!/usr/bin/env nextflow

process runColorSV {

    label 'colorsv'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(gfa_file), val(tumor_ids)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(filter_bed_file)
        val(params_colorsv_preprocess)
        val(params_colorsv_call)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_colorsv_outputs/"), emit: f

    script:
        """
        colorSV preprocess \
            -o ${sample_id}_colorsv_outputs \
            --graph $gfa_file \
            --reference $reference_genome_fasta_file \
            --tumor-ids $tumor_ids \
            $params_colorsv_preprocess

        colorSV call \
            -o ${sample_id}_colorsv_outputs \
            --graph $gfa_file \
            --filter $filter_bed_file \
            $params_colorsv_call
        """
}
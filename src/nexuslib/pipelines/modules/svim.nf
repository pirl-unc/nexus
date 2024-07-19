#!/usr/bin/env nextflow

process runSvimAlignmentMode {

    label 'svim'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(reference_genome_fasta_file)
        val(params_svim)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_svim.vcf"), emit: f

    script:
        """
        mkdir ${sample_id}_svim_output/
        svim alignment \
            --sample $sample_id \
            $params_svim \
            ${sample_id}_svim_output/ \
            $bam_file \
            $reference_genome_fasta_file
        cp ${sample_id}_svim_output/variants.vcf ${sample_id}_svim.vcf
        """
}
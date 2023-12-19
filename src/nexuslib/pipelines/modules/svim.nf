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
        val(svim)
        val(svim_params)
        val(output_dir)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_svim.vcf"), emit: f

    script:
        def svim_params_ = svim_params == true ? '' : svim_params

        """
        mkdir -p output/

        $svim alignment \
            --sample $sample_id \
            $svim_params_ \
            output/ \
            $bam_file \
            $reference_genome_fasta_file

        cp output/variants.vcf ${bam_file.baseName}_svim.vcf
        """
}
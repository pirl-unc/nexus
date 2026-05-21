#!/usr/bin/env nextflow

process runSvimAlignmentMode {

    label 'svim'
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
        val(params_svim)
        val(output_dir)

    output:
        tuple val(sample_id), path("svim/"), emit: f

    script:
        """
        mkdir svim/
        svim alignment \
            --sample $sample_id \
            $params_svim \
            svim/ \
            $bam_file \
            $reference_genome_fasta_file
        """
}
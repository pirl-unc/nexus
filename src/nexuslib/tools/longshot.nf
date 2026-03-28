#!/usr/bin/env nextflow

process runLongshot {

    label 'longshot'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_longshot)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_longshot.vcf"), emit: f

    script:
        """
        longshot \
            --bam $bam_file \
            --ref $reference_genome_fasta_file \
            --out ${sample_id}_longshot.vcf \
            $params_longshot
        """
}

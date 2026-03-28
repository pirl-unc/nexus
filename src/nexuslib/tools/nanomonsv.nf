#!/usr/bin/env nextflow

process runNanomonsv {

    label 'nanomonsv'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_nanomonsv_parse)
        val(params_nanomonsv_get)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_nanomonsv_outputs/"), emit: f

    script:
        """
        mkdir ${sample_id}_nanomonsv_outputs/
        nanomonsv parse \
            --reference_fasta $reference_genome_fasta_file \
            $params_nanomonsv_parse \
            $tumor_bam_file \
            ${sample_id}_nanomonsv_outputs/$sample_id
        nanomonsv parse \
            --reference_fasta $reference_genome_fasta_file \
            $params_nanomonsv_parse \
            $normal_bam_file \
            ${sample_id}_nanomonsv_outputs/${normal_bam_file.baseName}_nanomonsv
        nanomonsv get \
            --control_prefix ${sample_id}_nanomonsv_outputs/${normal_bam_file.baseName}_nanomonsv \
            --control_bam $normal_bam_file \
            --processes ${task.cpus} \
            $params_nanomonsv_get \
            ${sample_id}_nanomonsv_outputs/$sample_id \
            $tumor_bam_file \
            $reference_genome_fasta_file
        """
}

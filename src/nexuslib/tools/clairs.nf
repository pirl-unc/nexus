#!/usr/bin/env nextflow

process runClairs {

    label 'clairs'
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
        val(params_clairs)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_clairs_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_clairs_outputs/
        run_clairs \
            --tumor_bam $tumor_bam_file \
            --normal_bam $normal_bam_file \
            --ref_fn $reference_genome_fasta_file \
            --output_dir ${sample_id}_clairs_outputs/ \
            --threads ${task.cpus} \
            $params_clairs
        """
}
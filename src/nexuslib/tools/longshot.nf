#!/usr/bin/env nextflow

process runLongshot {

    label 'longshot'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy',
        pattern: "${sample_id}_longshot.vcf"
    )

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy',
        pattern: "${sample_id}_longshot.bam"
    )

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy',
        pattern: "${sample_id}_longshot.bam.bai"
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_longshot)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_longshot.vcf"), path("${sample_id}_longshot.bam"), path("${sample_id}_longshot.bam.bai"), emit: f

    script:
        """
        longshot \
            --bam $bam_file \
            --ref $reference_genome_fasta_file \
            --out ${sample_id}_longshot.vcf \
            --out_bam ${sample_id}_longshot.bam \
            $params_longshot
        samtools index -@ ${task.cpus} -b ${sample_id}_longshot.bam ${sample_id}_longshot.bam.bai
        """
}

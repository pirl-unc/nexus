#!/usr/bin/env nextflow

process runGridss {

    label 'gridss'
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
        val(params_gridss)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_gridss.vcf"), emit: f

    script:
        """
        gridss \
            --reference $reference_genome_fasta_file \
            --output ${sample_id}_gridss.vcf \
            --threads ${task.cpus} \
            --jvmheap ${task.memory.toGiga()}g \
            --otherjvmheap ${task.memory.toGiga()}g \
            $params_gridss \
            $normal_bam_file \
            $tumor_bam_file
        """
}

process runGridssSomaticFilter {

    label 'gridss'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(vcf_file)
        val(params_gridss_somatic_filter)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_gridss_high_confidence_somatic.vcf.bgz"), path("${sample_id}_gridss_high_confidence_somatic.vcf.bgz.tbi"), path("${sample_id}_gridss_high_and_low_confidence_somatic.vcf.bgz"), path("${sample_id}_gridss_high_and_low_confidence_somatic.vcf.bgz.tbi"), emit: f

    script:
        """
        Rscript /opt/gridss/gridss_somatic_filter \
            --input $vcf_file \
            --output ${sample_id}_gridss_high_confidence_somatic.vcf \
            --fulloutput ${sample_id}_gridss_high_and_low_confidence_somatic.vcf \
            --scriptdir /opt/gridss/ \
            -n 1 \
            -t 2 \
            $params_gridss_somatic_filter
        """
}
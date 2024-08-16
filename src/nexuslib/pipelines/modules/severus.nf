#!/usr/bin/env nextflow

process runSeverus {

    label 'severus'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), path(phased_vcf_file)
        path(vntr_bed_file)
        val(params_severus)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_severus_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_severus_outputs/
        severus \
            --target-bam $tumor_bam_file \
            --control-bam $normal_bam_file \
            --phasing-vcf $phased_vcf_file \
            --vntr-bed $vntr_bed_file \
            --threads ${task.cpus} \
            --out-dir ${sample_id}_severus_outputs/ \
            $params_severus
        """
}
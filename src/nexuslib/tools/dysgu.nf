#!/usr/bin/env nextflow

process runDysguSomaticMode {

    label 'dysgu'
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
        val(params_dysgu_run)
        val(params_dysgu_filter)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_dysgu_somatic_svs.vcf"), emit: f

    script:
        """
        dysgu run \
            --out-format vcf \
            --procs ${task.cpus} \
            $params_dysgu_run \
            $reference_genome_fasta_file \
            ${sample_id}_dysgu_run_tumor_temp/ \
            $tumor_bam_file > ${tumor_bam_file.baseName}.vcf
        dysgu run \
            --out-format vcf \
            --procs ${task.cpus} \
            $params_dysgu_run \
            $reference_genome_fasta_file \
            ${sample_id}_dysgu_run_normal_temp/ \
            $normal_bam_file > ${normal_bam_file.baseName}.vcf
        dysgu filter \
            --normal-vcf ${normal_bam_file.baseName}.vcf \
            --procs ${task.cpus} \
            $params_dysgu_filter \
            ${tumor_bam_file.baseName}.vcf \
            $normal_bam_file > ${sample_id}_dysgu_somatic_svs.vcf
        """
}

process runDysguGermlineMode {

    label 'dysgu'
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
        val(params_dysgu_run)
        val(params_dysgu_filter)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_dysgu_germline_svs.vcf"), emit: f

    script:
        """
        dysgu run \
            --out-format vcf \
            --procs ${task.cpus} \
            $params_dysgu_run \
            $reference_genome_fasta_file \
            ${sample_id}_dysgu_run_temp/ \
            $bam_file > ${bam_file.baseName}.vcf
        dysgu filter \
            --procs ${task.cpus} \
            $params_dysgu_filter \
            ${bam_file.baseName}.vcf > ${sample_id}_dysgu_germline_svs.vcf
        """
}

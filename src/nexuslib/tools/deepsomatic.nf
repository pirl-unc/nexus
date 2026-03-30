#!/usr/bin/env nextflow

process runDeepSomaticSingularity {

    label 'deepsomatic'
    tag "${sample_id}"
    scratch true
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        val(reference_genome_fasta_file)
        val(reference_genome_fasta_fai_file)
        val(singularity)
        val(deepsomatic_model_type)
        val(deepsomatic_bin_version)
        val(deepsomatic_bin_path)
        val(deepsomatic_input_path)
        val(deepsomatic_output_path)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_deepsomatic.vcf.gz"), path("${sample_id}_deepsomatic.gvcf.gz"), emit: f

    script:
        """
        mkdir -p \${PWD}/${sample_id}_deepsomatic_intermediate_results/
        $singularity run \
            -B ${deepsomatic_input_path}:${deepsomatic_input_path} \
            -B ${deepsomatic_output_path}:${deepsomatic_output_path} \
            docker://google/deepsomatic:$deepsomatic_bin_version \
            $deepsomatic_bin_path \
            --model_type=$deepsomatic_model_type \
            --ref=$reference_genome_fasta_file \
            --reads_tumor=$tumor_bam_file \
            --reads_normal=$normal_bam_file \
            --sample_name_tumor="tumor" \
            --sample_name_normal="normal" \
            --output_vcf=${sample_id}_deepsomatic.vcf.gz \
            --output_gvcf=${sample_id}_deepsomatic.gvcf.gz \
            --intermediate_results_dir=\${PWD}/${sample_id}_deepsomatic_intermediate_results/ \
            --num_shards=${task.cpus}
        """
}

process runDeepSomaticDocker {

    label 'deepsomatic'
    tag "${sample_id}"
    scratch true
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        val(reference_genome_fasta_file)
        val(reference_genome_fasta_fai_file)
        val(deepsomatic_model_type)
        val(deepsomatic_bin_version)
        val(deepsomatic_bin_path)
        val(deepsomatic_input_path)
        val(deepsomatic_output_path)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_deepsomatic.vcf.gz"), emit: f

    script:
        """
        mkdir -p \${PWD}/${sample_id}_deepsomatic_intermediate_results/
        docker run \
            -v ${deepsomatic_input_path}:${deepsomatic_input_path} \
            -v ${deepsomatic_output_path}:${deepsomatic_output_path} \
            google/deepsomatic:${deepsomatic_bin_version} \
            $deepsomatic_bin_path \
            --model_type=${deepsomatic_model_type} \
            --ref=$reference_genome_fasta_file \
            --reads_tumor=\${PWD}/$tumor_bam_file \
            --reads_normal=\${PWD}/$normal_bam_file \
            --sample_name_tumor="tumor" \
            --sample_name_normal="normal" \
            --output_vcf=\${PWD}/${sample_id}_deepsomatic.vcf.gz \
            --intermediate_results_dir=\${PWD}/${sample_id}_deepsomatic_intermediate_results/ \
            --num_shards=${task.cpus}
        """
}
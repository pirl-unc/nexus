#!/usr/bin/env nextflow

process runDeepVariantSingularity {

    label 'deepvariant'
    tag "${sample_id}"
    scratch true
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(reference_genome_fasta_file)
        val(singularity)
        val(deepvariant_model_type)
        val(deepvariant_bin_version)
        val(deepvariant_bin_path)
        val(deepvariant_input_path)
        val(deepvariant_output_path)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_deepvariant.vcf"), path("${sample_id}_deepvariant.gvcf"), emit: f

    script:
        """
        mkdir -p \${PWD}/${sample_id}_deepvariant_intermediate_results/
        $singularity run \
            -B ${deepvariant_input_path}:${deepvariant_input_path} \
            -B ${deepvariant_output_path}:${deepvariant_output_path} \
            docker://google/deepvariant:$deepvariant_bin_version \
            $deepvariant_bin_path \
            --model_type=$deepvariant_model_type \
            --ref=$reference_genome_fasta_file \
            --reads=$bam_file \
            --output_vcf=${sample_id}_deepvariant.vcf \
            --output_gvcf=${sample_id}_deepvariant.gvcf \
            --intermediate_results_dir=\${PWD}/${sample_id}_deepvariant_intermediate_results/ \
            --num_shards=${task.cpus}
        """
}

process runDeepVariantDocker {

    label 'deepvariant'
    tag "${sample_id}"
    scratch true
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), val(bam_file), val(bam_bai_file)
        val(reference_genome_fasta_file)
        val(deepvariant_model_type)
        val(deepvariant_bin_version)
        val(deepvariant_bin_path)
        val(deepvariant_input_path)
        val(deepvariant_output_path)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_deepvariant.vcf"), path("${sample_id}_deepvariant.gvcf"), emit: f

    script:
        """
        mkdir -p \${PWD}/${sample_id}_deepvariant_intermediate_results/
        docker run \
            -v ${deepvariant_input_path}:${deepvariant_input_path} \
            -v ${deepvariant_output_path}:${deepvariant_output_path} \
            google/deepvariant:${deepvariant_bin_version} \
            $deepvariant_bin_path \
            --model_type=${deepvariant_model_type} \
            --ref=$reference_genome_fasta_file \
            --reads=$bam_file \
            --output_vcf=\${PWD}/${sample_id}_deepvariant.vcf \
            --output_gvcf=\${PWD}/${sample_id}_deepvariant.gvcf \
            --intermediate_results_dir=\${PWD}/${sample_id}_deepvariant_intermediate_results/ \
            --num_shards=${task.cpus}
        """
}
#!/usr/bin/env nextflow

process runSVisionPro {

    label 'svisionpro'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        val(reference_genome_fasta_file)
        val(svisionpro_params)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}*.vcf"), path("${sample_id}*.log"), emit: f

    script:
        """
        mkdir -p output/
        SVision-pro \
            --target_path $tumor_bam_file \
            --base_path $normal_bam_file \
            --genome_path $reference_genome_fasta_file \
            --model_path /SVision-pro-1.8/src/pre_process/model_liteunet_1024_8_16_32_32_32.pth \
            --out_path output/ \
            --sample_name $sample_id \
            --process_num ${task.cpus} \
            $svisionpro_params
        mv output/* .
        """
}
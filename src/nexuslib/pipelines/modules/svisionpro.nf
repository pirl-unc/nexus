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
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(svisionpro_model_file)
        val(params_svisionpro)
        val(params_svisionpro_extract)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}*.vcf"), path("${sample_id}*.log"), emit: f

    script:
        """
        mkdir -p ${sample_id}_svisionpro_outputs/
        SVision-pro \
            --target_path $tumor_bam_file \
            --base_path $normal_bam_file \
            --genome_path $reference_genome_fasta_file \
            --model_path $svisionpro_model_file \
            --out_path ${sample_id}_svisionpro_outputs/ \
            --sample_name $sample_id \
            --process_num ${task.cpus} \
            $params_svisionpro
        mv ${sample_id}_svisionpro_outputs/* .
        python /opt/SVision-pro/extract_op.py \
            --input_vcf ${sample_id}*.vcf \
            $params_svisionpro_extract
        """
}

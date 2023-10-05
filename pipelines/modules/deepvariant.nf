#!/usr/bin/env nextflow

process runDeepVariant {

    label 'deepvariant'
    tag "${sample_id}"
    scratch true
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(reference_genome_fasta_file)
        val(singularity)
        val(deepvariant_model_type)
        val(deepvariant_bin_version)
        val(deepvariant_bin_path)
        val(deepvariant_lib_path)
        val(output_dir)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_deepvariant.vcf"), path("${bam_file.baseName}_deepvariant.gvcf"), emit: f

    script:
        """
        $singularity run \
            -B ${deepvariant_lib_path} \
            docker://google/deepvariant:$deepvariant_bin_version \
            $deepvariant_bin_path \
            --model_type=$deepvariant_model_type \
            --ref=$reference_genome_fasta_file \
            --reads=$bam_file \
            --output_vcf=${bam_file.baseName}_deepvariant.vcf \
            --output_gvcf=${bam_file.baseName}_deepvariant.gvcf \
            --intermediate_results_dir=\${PWD}/${bam_file.baseName}_deepvariant_intermediate_results/ \
            --num_shards=${task.cpus}
        """
}
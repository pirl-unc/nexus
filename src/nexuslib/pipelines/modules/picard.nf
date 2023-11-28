#!/usr/bin/env nextflow

process runPicardMergeVCFs {

    label 'picard_mergevcfs'
    tag "${sample_id}"
    debug true

   publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(vcf_files)
        val(picard)
        val(variant_calling_method)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_${variant_calling_method}.vcf"), emit: f

    script:
        """
        vcf_files_str="${vcf_files}"
        IFS=' ' array=(\${vcf_files_str})
        vcf_files_params=""
        for element in "\${array[@]}"; do
            vcf_files_params+="I=\${element} "
        done
        java -Xmx64G -jar $picard MergeVcfs \
            \${vcf_files_params} \
            O=${sample_id}_${variant_calling_method}.vcf
        """
}
#!/usr/bin/env nextflow

process runPicardMergeVCFs {

    label 'picard_mergevcfs'
    tag "${sample_id}"
    debug true

   publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(vcf_files)
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
        java -Xmx${task.java_max_mem.toGiga()}G -jar /opt/picard/picard.jar MergeVcfs \
            \${vcf_files_params} \
            O=${sample_id}_${variant_calling_method}.vcf
        """
}

process runPicardSingleEndFastqToSam {

    label 'picard_fastqtosam'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_picard)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}.sam"), emit: f

    script:
        """
        java -Xmx${task.java_max_mem.toGiga()}G -jar /opt/picard/picard.jar FastqToSam \
            F1=${fastq_file} \
            SM=${sample_id} \
            O=${sample_id}.sam \
            $params_picard
        """
}

process runPicardPairedEndFastqToSam {

    label 'picard_fastqtosam'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(params_picard)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}.sam"), emit: f

    script:
        """
        java -Xmx${task.java_max_mem.toGiga()}G -jar /opt/picard/picard.jar FastqToSam \
            F1=${fastq_file_1} \
            F2=${fastq_file_2} \
            SM=${sample_id} \
            O=${sample_id}.sam \
            $params_picard
        """
}

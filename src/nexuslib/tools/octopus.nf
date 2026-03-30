#!/usr/bin/env nextflow

process runOctopusSomatic {

    label 'octopus'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(regions_txt_file)
        val(params_octopus)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_octopus.vcf"), emit: f

    script:
        """
        normal_sample_id=\$(samtools view -H "${normal_bam_file}" | awk -F'\t' '/^@RG/ { for(i=1;i<=NF;i++) if (\$i ~ /^SM:/) print substr(\$i,4) }')
        octopus \
            --reads $normal_bam_file $tumor_bam_file \
            --reference $reference_genome_fasta_file \
            --normal-sample "\${normal_sample_id}" \
            -o ${sample_id}_octopus.vcf \
            --threads ${task.cpus} \
            --target-working-memory ${task.per_thread_memory.toGiga()}G \
            --regions-file $regions_txt_file \
            $params_octopus
        """
}

process runOctopusGermline {

    label 'octopus'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(regions_txt_file)
        val(params_octopus)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_octopus.vcf"), emit: f

    script:
        """
        octopus \
            --reads $bam_file \
            --reference $reference_genome_fasta_file \
            --threads ${task.cpus} \
            --target-working-memory ${task.per_thread_memory.toGiga()}G \
            --regions-file $regions_txt_file \
            -o ${sample_id}_octopus.vcf \
            $params_octopus
        """
}

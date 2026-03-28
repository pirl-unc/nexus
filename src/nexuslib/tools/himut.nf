#!/usr/bin/env nextflow

process runHimutCall {

    label 'himut'
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
        path(region_list_file)
        val(params_himut)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_himut.vcf"), emit: f

    script:
        """
        himut call \
            --bam $bam_file \
            --ref $reference_genome_fasta_file \
            --region_list $region_list_file \
            --threads ${task.cpus} \
            --output ${sample_id}_himut.vcf \
            $params_himut
        """
}

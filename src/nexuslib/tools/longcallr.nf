#!/usr/bin/env nextflow

process runLongcallR {

    label 'longcallr'
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
        path(reference_genes_gtf_file)
        val(preset)
        val(params_longcallr)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_longcallr.phased.bam"), path("${sample_id}_longcallr.vcf"), emit: f

    script:
        """
        longcallR \
            --bam-path $bam_file \
            --ref-path $reference_genome_fasta_file \
            --annotation $reference_genes_gtf_file \
            --preset $preset \
            --output ${sample_id}_longcallr \
            --threads ${task.cpus} \
            $params_longcallr
        """
}
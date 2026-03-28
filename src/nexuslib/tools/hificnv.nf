#!/usr/bin/env nextflow

process runHiFiCNV {

    label 'hificnv'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(vcf_file)
        path(reference_genome_fasta_file)
        path(exclude_bed_file)
        val(params_hificnv)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_hificnv_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_hificnv_outputs/
        hificnv \
            --ref $reference_genome_fasta_file \
            --bam $bam_file \
            --maf $vcf_file \
            --exclude $exclude_bed_file \
            --output-prefix ${sample_id}_hificnv_outputs/${sample_id}_hificnv \
            --threads ${task.cpus} \
            $params_hificnv
        """
}
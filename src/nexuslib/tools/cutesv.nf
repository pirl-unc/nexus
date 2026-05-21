#!/usr/bin/env nextflow

process runCuteSV {

    label 'cutesv'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/cutesv/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_cutesv)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_cutesv.vcf"), emit: f

    script:
        """
        mkdir -p temp/
        cuteSV \
            --threads ${task.cpus} \
            --sample $sample_id \
            $params_cutesv \
            $bam_file \
            $reference_genome_fasta_file \
            ${sample_id}_cutesv.vcf \
            temp/
        """
}
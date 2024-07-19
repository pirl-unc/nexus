#!/usr/bin/env nextflow

process runSniffles2 {

    label 'sniffles2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(reference_genome_fasta_file)
        val(params_sniffles2)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_sniffles2.vcf"), emit: f

    script:
        """
        sniffles \
            --input $bam_file \
            --vcf ${sample_id}_sniffles2.vcf \
            --snf ${sample_id}_sniffles2.snf \
            --sample-id ${sample_id} \
            --reference ${reference_genome_fasta_file} \
            --threads ${task.cpus} \
            $params_sniffles2
        """
}

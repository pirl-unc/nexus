#!/usr/bin/env nextflow

process runMandalorion {

    label 'mandalorion'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genes_gtf_file)
        val(params_mandalorion)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_mandalorion_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_mandalorion_outputs/
        printf "%s\n" "$fastq_file" >> sample.fofn

        Mando.py \
            --path ${sample_id}_mandalorion_outputs/ \
            --genome_annotation $reference_genes_gtf_file \
            --genome_sequence $reference_genome_fasta_file \
            -f sample.fofn \
            -t ${task.cpus} \
            $params_mandalorion
        """
}

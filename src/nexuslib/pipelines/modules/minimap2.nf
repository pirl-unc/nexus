#!/usr/bin/env nextflow

process runMinimap2 {

    label 'minimap2'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_minimap2)
        val(platform_tag)
        val(platform_unit_tag)
        val(library_tag)

    output:
        tuple val(sample_id), path("${sample_id}_minimap2.sam"), emit: f

    script:
        """
        minimap2 \
            $params_minimap2 \
            -t ${task.cpus} \
            -R "@RG\\tID:$sample_id\\tSM:$sample_id\\tPL:$platform_tag\\tLB:$library_tag\\tPU:$platform_unit_tag" \
            $reference_genome_fasta_file $fastq_file > ${sample_id}_minimap2.sam
        """
}

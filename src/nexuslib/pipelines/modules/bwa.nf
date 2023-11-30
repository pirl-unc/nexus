#!/usr/bin/env nextflow

process runBwaMem2 {

    label 'bwa_mem2'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(bwa_mem2)
        val(samtools)
        val(reference_genome_fasta_file)
        val(platform_tag)
        val(platform_unit_tag)
        val(library_tag)

    output:
        tuple val(sample_id), path("${sample_id}.sam"), emit: f

    script:
        """
        $bwa_mem2 mem -t ${task.cpus} \
            -R "@RG\tID:$sample_id\tSM:$sample_id\tPL:$platform_tag\tLB:$library_tag\tPU:$platform_unit_tag" \
            $reference_genome_fasta_file \
            $fastq_file_1 \
            $fastq_file_2 > ${sample_id}.sam
        """
}

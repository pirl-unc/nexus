#!/usr/bin/env nextflow

process runUltraIndex {

    label 'ultra_index'
    debug true

    input:
        path(reference_genome_fasta_file)
        path(reference_gtf_file)
        val(params_ultra_index)

    output:
        path("ultra_index/"), emit: f

    script:
        """
        mkdir -p ultra_index
        uLTRA index \
            $params_ultra_index \
            $reference_genome_fasta_file \
            $reference_gtf_file \
            ultra_index/
        """
}

process runUltra {

    label 'ultra'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(ultra_index)
        val(params_ultra)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_ultra_sorted.bam"), path("${sample_id}_ultra_sorted.bam.bai"), emit: f

    script:
        """
        mkdir -p ${sample_id}_ultra/
        uLTRA align \
            --t ${task.cpus} \
            --index $ultra_index \
            --prefix ${sample_id}_ultra \
            $params_ultra \
            $reference_genome_fasta_file \
            $fastq_file \
            ${sample_id}_ultra/
        samtools sort -@ {task.cpus} -m ${task.samtools_memory.toGiga()}G -O bam -o ${sample_id}_ultra_sorted.bam ${sample_id}_ultra/minimap2.sam
        samtools index -b ${sample_id}_ultra_sorted.bam ${sample_id}_ultra_sorted.bam.bai
        """
}
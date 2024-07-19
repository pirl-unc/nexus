#!/usr/bin/env nextflow

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
        gunzip -c $fastq_file > ${sample_id}.fastq
        uLTRA align \
            --t ${task.cpus} \
            --index $ultra_index \
            --prefix ${sample_id}_ultra \
            $params_ultra \
            $reference_genome_fasta_file \
            ${sample_id}.fastq \
            ${sample_id}_ultra/
        samtools sort -@ {task.cpus} -m ${task.samtools_memory.toGiga()}G -O bam -o ${sample_id}_ultra_sorted.bam ${sample_id}_ultra/minimap2.sam
        samtools index -b ${sample_id}_ultra_sorted.bam ${sample_id}_ultra_sorted.bam.bai
        """
}
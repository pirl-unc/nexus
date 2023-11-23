#!/usr/bin/env nextflow

process runUltra {

    label 'ultra'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file)
        val(reference_genome_fasta_file)
        val(ultra)
        val(ultra_index)
        val(ultra_params)
        val(samtools)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_ultra_sorted.bam"), path("${sample_id}_ultra_sorted.bam.bai"), emit: f

    script:
        """
        mkdir -p ${sample_id}_ultra/
        gunzip -c $fastq_file > ${sample_id}.fastq
        $ultra align \
            --t ${task.ultra_cpus} \
            --index $ultra_index \
            --prefix ${sample_id}_ultra \
            $ultra_params \
            $reference_genome_fasta_file \
            ${sample_id}.fastq \
            ${sample_id}_ultra/
        $samtools sort -@ {task.samtools_cpus} -m ${task.samtools_memory.toGiga()}G -O bam -o ${sample_id}_ultra_sorted.bam ${sample_id}_ultra/minimap2.sam
        $samtools index -b ${sample_id}_ultra_sorted.bam ${sample_id}_ultra_sorted.bam.bai
        """
}
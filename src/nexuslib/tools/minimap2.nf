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

process runMinimap2SortedBam {

    label 'minimap2_sorted_bam'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(calmd_reference_fasta_file)
        path(calmd_reference_fasta_fai_file)
        val(params_minimap2)
        val(platform_tag)
        val(platform_unit_tag)
        val(library_tag)

    output:
        tuple val(sample_id), path("${sample_id}_minimap2_sorted.bam"), path("${sample_id}_minimap2_sorted.bam.bai"), emit: f

    script:
        """
        minimap2 \
            $params_minimap2 \
            -t ${task.minimap2_threads} \
            -R "@RG\\tID:$sample_id\\tSM:$sample_id\\tPL:$platform_tag\\tLB:$library_tag\\tPU:$platform_unit_tag" \
            $reference_genome_fasta_file $fastq_file \
            | samtools view -@ ${task.samtools_view_threads} -bS \
            | samtools calmd -@ ${task.samtools_calmd_threads} -b - $calmd_reference_fasta_file \
            | samtools sort -@ ${task.samtools_sort_threads} -m ${task.samtools_memory.toGiga()}G -O bam -o ${sample_id}_minimap2_sorted.bam
        samtools index -@ ${task.cpus} -b ${sample_id}_minimap2_sorted.bam ${sample_id}_minimap2_sorted.bam.bai
        """
}

process runMinimap2CustomReference {

    label 'minimap2'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file), path(reference_genome_fasta_file), path(reference_genome_fasta_fai_file)
        val(params_minimap2)
        val(platform_tag)
        val(platform_unit_tag)
        val(library_tag)

    output:
        tuple val(sample_id), path("${sample_id}_minimap2.sam"), path(reference_genome_fasta_file), path(reference_genome_fasta_fai_file), emit: f

    script:
        """
        minimap2 \
            $params_minimap2 \
            -t ${task.cpus} \
            -R "@RG\\tID:$sample_id\\tSM:$sample_id\\tPL:$platform_tag\\tLB:$library_tag\\tPU:$platform_unit_tag" \
            $reference_genome_fasta_file $fastq_file > ${sample_id}_minimap2.sam
        """
}

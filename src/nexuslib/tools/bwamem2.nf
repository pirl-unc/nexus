#!/usr/bin/env nextflow

process runBwaMem2Index {

    label 'bwamem2_index'
    debug true

    input:
        path(fasta_file)

    output:
        path("bwamem2_index"), emit: f

    script:
        """
        mkdir -p bwamem2_index
        bwa-mem2 index -p bwamem2_index/$fasta_file $fasta_file
        """
}

process runBwaMem2 {

    label 'bwamem2'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_dict_file)
        path(reference_genome_fasta_0123_file)
        path(reference_genome_fasta_amb_file)
        path(reference_genome_fasta_ann_file)
        path(reference_genome_fasta_bwt_file)
        path(reference_genome_fasta_pac_file)
        val(platform_tag)
        val(platform_unit_tag)
        val(library_tag)

    output:
        tuple val(sample_id), path("${sample_id}_bwamem2_sorted.bam"), path("${sample_id}_bwamem2_sorted.bam.bai"), emit: f

    script:
        """
        bwa-mem2 mem -t ${task.bwamem2_threads} \\
            -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:${platform_tag}\\tLB:${library_tag}\\tPU:${platform_unit_tag}" \\
            ${reference_genome_fasta_file} ${fastq_file_1} ${fastq_file_2} \\
            | samtools view -@ ${task.samtools_view_threads} -bS \\
            | samtools sort -@ ${task.samtools_sort_threads} -m ${task.samtools_memory.toGiga()}G -O bam -o ${sample_id}_bwamem2_sorted.bam
        samtools index -@ ${task.cpus} -b ${sample_id}_bwamem2_sorted.bam ${sample_id}_bwamem2_sorted.bam.bai
        """
}

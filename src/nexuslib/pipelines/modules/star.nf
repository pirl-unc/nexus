#!/usr/bin/env nextflow

process runStar {

    label 'star'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        path(star_index)
        val(params_star)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_star_*.bam"), path("${sample_id}_star_*.bam.bai"), emit: f

    script:
        """
        STAR \
            --runThreadN ${task.cpus} \
            --readFilesIn $fastq_file_1 $fastq_file_2 \
            --genomeDir $star_index \
            --outFileNamePrefix ${sample_id}_star_ \
            $params_star
        samtools index -b ${sample_id}_star_*.bam
        """
}

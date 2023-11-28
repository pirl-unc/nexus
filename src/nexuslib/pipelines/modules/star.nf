#!/usr/bin/env nextflow

process runStar {

    label 'star'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(star)
        val(star_params)
        val(star_index)

    output:
        tuple val(sample_id), path("${sample_id}_star_*.bam"), emit: f

    script:
        """
        $star \
            --runThreadN ${task.cpus} \
            --readFilesIn $fastq_file_1 $fastq_file_2 \
            --genomeDir $star_index \
            --outFileNamePrefix ${sample_id}_star_ \
            $star_params
        """
}

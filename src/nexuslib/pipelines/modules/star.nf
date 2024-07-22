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
        tuple val(sample_id), emit: f

    script:
        """
        STAR \
            --runThreadN ${task.cpus} \
            --readFilesIn $fastq_file_1 $fastq_file_2 \
            --genomeDir $star_index \
            $params_star
        """
}

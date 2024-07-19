#!/usr/bin/env nextflow

process runRatatoskCorrectFirstPass {

    label 'ratatosk'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(long_read_fastq_file), path(short_read_r1_fastq_file), path(short_read_r2_fastq_file)
        val(params_ratatosk)

    output:
        tuple val(sample_id), path(long_read_fastq_file), path(short_read_r1_fastq_file), path(short_read_r2_fastq_file), path("${long_read_fastq_file.simpleName}.2.fastq"), emit: f

    script:
        """
        Ratatosk correct \
            -1 \
            -v \
            -c ${task.cpus} \
            -l $long_read_fastq_file \
            -s $short_read_r1_fastq_file \
            -s $short_read_r2_fastq_file \
            -o ${long_read_fastq_file.simpleName} \
            $params_ratatosk
        """
}

process runRatatoskCorrectSecondPass {

    label 'ratatosk'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(long_read_fastq_file), path(short_read_r1_fastq_file), path(short_read_r2_fastq_file), path(long_read_corrected_fastq_file)
        val(params_ratatosk)
        val(output_dir)

    output:
        tuple val(sample_id), path("${long_read_fastq_file.simpleName}_corrected.fastq.gz"), emit: f

    script:
        """
        Ratatosk correct \
            -2 \
            -v \
            -G \
            -c ${task.cpus} \
            -l $long_read_corrected_fastq_file \
            -s $short_read_r1_fastq_file \
            -s $short_read_r2_fastq_file \
            -L $long_read_fastq_file \
            -o ${long_read_fastq_file.simpleName}_corrected \
            $params_ratatosk
        """
}

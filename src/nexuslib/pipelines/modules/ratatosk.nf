#!/usr/bin/env nextflow

process runRatatoskIndexFirstPass {

    label 'ratatosk'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(long_read_fastq_file), path(short_read_r1_fastq_file), path(short_read_r2_fastq_file)
        val(ratatosk)
        val(ratatosk_params)

    output:
        tuple val(sample_id), path(long_read_fastq_file), path("${long_read_fastq_file.simpleName}_pass1.index.k31.fasta"), path("${long_read_fastq_file.simpleName}_pass1.index.k31.rtsk"), path("${long_read_fastq_file.simpleName}_pass1.index.k63.fasta"), emit: f

    script:
        def ratatosk_params_ = ratatosk_params == true ? '' : ratatosk_params

        """
        $ratatosk index \
            -1 \
            -v \
            -c ${task.cpus} \
            -l $long_read_fastq_file \
            -s $short_read_r1_fastq_file \
            -s $short_read_r2_fastq_file \
            -o ${long_read_fastq_file.simpleName}_pass1 \
            $ratatosk_params_
        """
}

process runRatatoskCorrectFirstPass {

    label 'ratatosk'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(long_read_fastq_file), path(k31_fasta_file), path(k31_rtsk_file), path(k63_fasta_file)
        val(ratatosk)
        val(ratatosk_params)

    output:
        tuple val(sample_id), path("${long_read_fastq_file.simpleName}_pass1.2.fastq"), path(long_read_fastq_file), path(k63_fasta_file), emit: f

    script:
        def ratatosk_params_ = ratatosk_params == true ? '' : ratatosk_params

        """
        $ratatosk correct \
            -1 \
            -v \
            -c ${task.cpus} \
            -g ${k31_fasta_file} \
            -d ${k31_rtsk_file} \
            -l $long_read_fastq_file \
            -o ${long_read_fastq_file.simpleName}_pass1 \
            $ratatosk_params_
        """
}

process runRatatoskIndexSecondPass {

    label 'ratatosk'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(long_read_pass1_corrected_fastq_file), path(long_read_original_fastq_file), path(k63_fasta_file)
        val(ratatosk)
        val(ratatosk_params)

    output:
        tuple val(sample_id), path(long_read_pass1_corrected_fastq_file), path(long_read_original_fastq_file), path("${long_read_original_fastq_file.simpleName}_pass2.index.k63.fasta"), path("${long_read_original_fastq_file.simpleName}_pass2.index.k63.rtsk"), emit: f

    script:
        def ratatosk_params_ = ratatosk_params == true ? '' : ratatosk_params

        """
        $ratatosk index \
            -2 \
            -v \
            -c ${task.cpus} \
            -g $k63_fasta_file \
            -l $long_read_pass1_corrected_fastq_file \
            -o ${long_read_original_fastq_file.simpleName}_pass2 \
            $ratatosk_params_
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
        tuple val(sample_id), path(long_read_pass1_corrected_fastq_file), path(long_read_original_fastq_file), path(k63_fasta_file), path(k63_rtsk_file)
        val(ratatosk)
        val(ratatosk_params)
        val(output_dir)

    output:
        tuple val(sample_id), path("${long_read_original_fastq_file.simpleName}_ratatosk_corrected.fastq.gz"), emit: f

    script:
        def ratatosk_params_ = ratatosk_params == true ? '' : ratatosk_params

        """
        $ratatosk correct \
            -2 \
            -v \
            -c ${task.cpus} \
            -g ${k63_fasta_file} \
            -d ${k63_rtsk_file} \
            -l ${long_read_pass1_corrected_fastq_file} \
            -L ${long_read_original_fastq_file} \
            -o ${long_read_original_fastq_file.simpleName}_ratatosk_corrected \
            $ratatosk_params_
        gzip ${long_read_original_fastq_file.simpleName}_ratatosk_corrected.fastq
        """
}

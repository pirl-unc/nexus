#!/usr/bin/env nextflow

process runSalmonAlignmentMode {

    label 'salmon'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(salmon)
        val(salmon_params)
        val(output_dir)

    output:
        tuple val(sample_id), path("salmon/"), emit: f

    script:
        def salmon_params_ = salmon_params == true ? '' : salmon_params

        """
        mkdir -p salmon/
        $salmon quant \
            --alignments $bam_file \
            --output salmon/ \
            --threads ${task.cpus} \
            $salmon_params_
        """
}

process runSalmonPairedEndMappingMode {

    label 'salmon'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2), path(transcripts_fasta_file)
        val(salmon)
        val(salmon_index_params)
        val(salmon_quant_params)
        val(output_dir)

    output:
        tuple val(sample_id), path("salmon/"), emit: f

    script:
        def salmon_index_params_ = salmon_index_params == true ? '' : salmon_index_params
        def salmon_quant_params_ = salmon_quant_params == true ? '' : salmon_quant_params

        """
        mkdir -p index/
        $salmon index \
            $salmon_index_params_ \
            --transcripts $transcripts_fasta_file \
            --index index/ \
            --threads ${task.cpus}

        mkdir -p salmon/
        $salmon quant \
            $salmon_quant_params_ \
            --index index/ \
            --mates1 $fastq_file_1 \
            --mates2 $fastq_file_2 \
            --output salmon/ \
            --threads ${task.cpus}
        """
}


#!/usr/bin/env nextflow

process runFlair3Transcriptome {

    label 'flair3'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/flair3/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genes_gtf_file)
        val(params_flair)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_flair3.*"), emit: f

    script:
        """
        flair transcriptome \
          -b $bam_file \
          --genome $reference_genome_fasta_file \
          --gtf $reference_genes_gtf_file \
          --output ${sample_id}_flair3 \
          --threads ${task.cpus} \
          $params_flair
        """
}

process runFlair3Align {

    label 'flair3'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/flair3/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        path(reference_genome_fasta_file)
        val(params_flair)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_flair3.bam"), path("${sample_id}_flair3.bam.bai"), path("${sample_id}_flair3.bed"), emit: f

    script:
        """
        flair align \
          -g $reference_genome_fasta_file \
          -r $fastq_file \
          --output ${sample_id}_flair3 \
          --threads ${task.cpus} \
          $params_flair
        """
}
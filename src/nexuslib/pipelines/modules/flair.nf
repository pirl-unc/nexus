#!/usr/bin/env nextflow

process runFlairAlign {

    label 'flair'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        path(reference_genome_fasta_file)
        val(params_flair)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), path("${sample_id}.bed"), emit: f

    script:
        """
        flair align \
          -g $reference_genome_fasta_file \
          -r $fastq_file \
          --output $sample_id \
          --threads ${task.cpus} \
          $params_flair
        """
}

process runFlairCorrect {

    label 'flair'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(bed_file)
        path(reference_genome_fasta_file)
        path(gtf_file)
        val(params_flair)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_all_corrected.bed"), path("${sample_id}_all_inconsistent.bed"), emit: f

    script:
        """
        flair correct \
          --query $bed_file \
          --genome $reference_genome_fasta_file \
          --gtf $gtf_file \
          --output $sample_id \
          --threads ${task.cpus} \
          $params_flair
        """
}

process runFlairCollapse {

    label 'flair'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(all_corrected_bed), path(all_inconsistent_bed), path(fastq_file)
        path(reference_genome_fasta_file)
        path(gtf_file)
        val(params_flair)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}.isoforms.bed"), path("${sample_id}.isoforms.fa"), path("${sample_id}.isoforms.gtf"), emit: f

    script:
        """
        flair collapse \
          --query $all_corrected_bed \
          --genome $reference_genome_fasta_file \
          --gtf $gtf_file \
          --reads $fastq_file \
          --output $sample_id \
          --threads ${task.cpus} \
          $params_flair
        """
}

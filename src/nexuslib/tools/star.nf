#!/usr/bin/env nextflow

process runStarIndex {

    label 'star_index'
    debug true

    input:
        path(reference_genome_fasta_file)
        path(reference_gtf_file)
        val(params_star_genomegenerate)

    output:
        path("star_index/"), emit: f

  script:
      """
      mkdir -p star_index/
      STAR \
          --runMode genomeGenerate \
          --genomeDir star_index/ \
          --genomeFastaFiles $reference_genome_fasta_file \
          --sjdbGTFfile $reference_gtf_file \
          --runThreadN ${task.cpus} \
          $params_star_genomegenerate
      """
}

process runStar {

    label 'star'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        path(star_index)
        val(params_star)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_star_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_star_outputs/
        STAR \
            --runThreadN ${task.cpus} \
            --readFilesIn $fastq_file_1 $fastq_file_2 \
            --genomeDir $star_index \
            --outFileNamePrefix ${sample_id}_star_outputs/${sample_id}_star_ \
            $params_star
        bam_file=\$(find ${sample_id}_star_outputs/ -name "${sample_id}_star_*.bam" | head -n 1)
        samtools index -b \$bam_file
        """
}

process runStarNoPublish {

    label 'star'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        path(star_index)
        val(params_star)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_star_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_star_outputs/
        STAR \
            --runThreadN ${task.cpus} \
            --readFilesIn $fastq_file_1 $fastq_file_2 \
            --genomeDir $star_index \
            --outFileNamePrefix ${sample_id}_star_outputs/${sample_id}_star_ \
            $params_star
        bam_file=\$(find ${sample_id}_star_outputs/ -name "${sample_id}_star_*.bam" | head -n 1)
        samtools index -b \$bam_file
        """
}

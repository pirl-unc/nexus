#!/usr/bin/env nextflow

process runIsoquant {

    label 'isoquant'
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        path(reference_genome_fasta_file)
        path(gtf_file)
        val(params_isoquant)
        val(output_dir)

    output:
        tuple val(sample_id), path("$sample_id/"), emit: f

    script:
        """
        echo "#${sample_id}" > fastq_list.txt
        echo "$fastq_file:$sample_id" >> fastq_list.txt
        output_dirname=$sample_id
        mkdir -p $sample_id/
        isoquant.py \
          --reference $reference_genome_fasta_file \
          --genedb $gtf_file \
          --fastq_list fastq_list.txt \
          --threads ${task.cpus} \
          -o $sample_id/ \
          $params_isoquant
        """
}


#!/usr/bin/env nextflow

process runDiamondBlastp {

    label 'diamond_blastp'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fasta_file)
        path(database_file)
        val(params_diamond_blastp)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_diamond_blastp_output.tsv"), emit: f

    script:
        """
        diamond blastp \
            --query $fasta_file \
            --db $database_file \
            --out ${sample_id}_diamond_blastp_output.tsv \
            --header verbose \
            --verbose \
            --threads ${task.cpus} \
            $params_diamond_blastp
        """
}

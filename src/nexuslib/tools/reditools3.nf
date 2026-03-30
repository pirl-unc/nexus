#!/usr/bin/env nextflow

process runReditools3Analyze {

    label 'reditools3'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_reditools3_analyze)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_reditools3_analyze.tsv"), emit: f

    script:
        """
        python -m reditools analyze \
            --threads ${task.cpus} \
            --reference $reference_genome_fasta_file \
            --output-file ${sample_id}_reditools3_analyze.tsv \
            $params_reditools3_analyze \
            $bam_file
        """
}
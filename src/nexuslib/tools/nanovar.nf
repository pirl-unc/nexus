#!/usr/bin/env nextflow

process runNanoVar {

    label 'nanovar'
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
        val(params_nanovar)
        val(output_dir)

    output:
        tuple val(sample_id), path("nanovar/"), emit: f

    script:
        """
        mkdir -p nanovar/
        nanovar \
            $params_nanovar \
            --threads ${task.cpus} \
            $bam_file \
            $reference_genome_fasta_file \
            nanovar/
        """
}

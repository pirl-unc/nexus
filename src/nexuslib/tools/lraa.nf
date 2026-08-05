#!/usr/bin/env nextflow

process runLRAA {

    label 'lraa'
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
        path(gtf_file)
        val(params_lraa)
        val(output_dir)

    output:
        tuple val(sample_id), path("lraa/"), emit: f

    script:
        """
        mkdir lraa/
        LRAA \
            --genome $reference_genome_fasta_file \
            --gtf $gtf_file \
            --bam $bam_file \
            --output_prefix lraa/${sample_id} \
            --num_threads_per_worker ${Math.max(1, task.cpus.intdiv(4))} \
            $params_lraa
        """
}

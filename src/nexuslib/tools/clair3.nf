#!/usr/bin/env nextflow

process runClair3 {

    label 'clair3'
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
        val(params_clair3)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_clair3_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_clair3_outputs/
        run_clair3.sh \
            --bam_fn=$bam_file \
            --ref_fn=$reference_genome_fasta_file \
            --output=${sample_id}_clair3_outputs/ \
            --threads ${task.cpus} \
            $params_clair3
        """
}

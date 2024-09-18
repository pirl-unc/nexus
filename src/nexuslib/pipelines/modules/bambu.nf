#!/usr/bin/env nextflow

process runBambuDiscoverQuantify {

    label 'bambu'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(gtf_file)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_bambu_discover_quantify_output/"), emit: f

    script:
        """
        mkdir ${sample_id}_bambu_discover_quantify_output/
        Rscript /opt/bambu/run_bambu.R \
            --bam-file $bam_file \
            --fasta-file $reference_genome_fasta_file \
            --gtf-file $gtf_file \
            --output-path ${sample_id}_bambu_discover_quantify_output/
        """
}

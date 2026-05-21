#!/usr/bin/env nextflow

process runSvaba {

    label 'svaba'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_gzi_file)
        path(reference_genome_fasta_amb_file)
        path(reference_genome_fasta_ann_file)
        path(reference_genome_fasta_bwt_file)
        path(reference_genome_fasta_pac_file)
        path(reference_genome_fasta_sa_file)
        val(params_svaba)
        val(output_dir)

    output:
        tuple val(sample_id), path("svaba/"), emit: f

    script:
        """
        mkdir -p svaba/
        svaba run \
            --reference-genome $reference_genome_fasta_file \
            --id-string svaba/$sample_id \
            --case-bam $tumor_bam_file \
            --control-bam $normal_bam_file \
            --threads ${task.cpus} \
            $params_svaba
        """
}
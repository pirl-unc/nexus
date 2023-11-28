#!/usr/bin/env nextflow

process runSavana {

    label 'savana'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        val(reference_genome_fasta_file)
        val(savana)
        val(savana_params)
        val(output_dir)

    output:
        tuple val(sample_id), path("${tumor_bam_file.baseName}.sv_breakpoints.vcf"), path("${tumor_bam_file.baseName}.sv_breakpoints_read_support.tsv"), path("${tumor_bam_file.baseName}.sv_breakpoints.bedpe"), emit: f

    script:
        def savana_params_ = savana_params == true ? '' : savana_params

        """
        mkdir -p savana/
        $savana run \
            -t $tumor_bam_file \
            -n $normal_bam_file \
            --ref $reference_genome_fasta_file \
            --threads ${task.cpus} \
            --outdir savana/ \
            --sample ${tumor_bam_file.baseName} \
            $savana_params_
        cp savana/*.bedpe .
        cp savana/*.tsv .
        cp savana/*.vcf .
        """
}
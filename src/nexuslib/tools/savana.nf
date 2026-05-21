#!/usr/bin/env nextflow

process runSavanaRun {

    label 'savana'
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
        path(contigs_file)
        val(params_savana_run)
        val(output_dir)

    output:
        tuple val(sample_id), path("savana/"), emit: f

    script:
        """
        mkdir -p savana/
        savana run \
            -t $tumor_bam_file \
            -n $normal_bam_file \
            --ref $reference_genome_fasta_file \
            --ref_index $reference_genome_fasta_fai_file \
            --contigs $contigs_file \
            --threads ${task.cpus} \
            --outdir savana/ \
            --sample $sample_id \
            $params_savana_run
        """
}

process runSavanaClassify {

    label 'savana'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/savana/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(savana_run_dir)
        path(custom_params_file)
        val(params_savana_classify)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_savana_classify_output.somatic.vcf"), path("${sample_id}_savana_classify_output.germline.vcf"), emit: f

    script:
        """
        savana classify \
            --vcf ${savana_run_dir}/${sample_id}.sv_breakpoints.vcf \
            --output ${sample_id}_savana_classify_output.vcf \
            --custom_params $custom_params_file \
            $params_savana_classify
        """
}
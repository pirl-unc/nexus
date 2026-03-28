#!/usr/bin/env nextflow

process runSavanaRun {

    label 'savana'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
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
        tuple val(sample_id), path("${sample_id}_savana_run_breakpoints.vcf"), path("${sample_id}_savana_run_breakpoints_read_support.tsv"), path("${sample_id}_savana_run_breakpoints.bedpe"), emit: f

    script:
        """
        mkdir -p ${sample_id}_savana_run_outputs/
        savana run \
            -t $tumor_bam_file \
            -n $normal_bam_file \
            --ref $reference_genome_fasta_file \
            --ref_index $reference_genome_fasta_fai_file \
            --contigs $contigs_file \
            --threads ${task.cpus} \
            --outdir ${sample_id}_savana_run_outputs/ \
            --sample $sample_id \
            $params_savana_run
        mv ${sample_id}_savana_run_outputs/${sample_id}.sv_breakpoints.vcf ${sample_id}_savana_run_breakpoints.vcf
        mv ${sample_id}_savana_run_outputs/${sample_id}.sv_breakpoints.bedpe ${sample_id}_savana_run_breakpoints.bedpe
        mv ${sample_id}_savana_run_outputs/${sample_id}.sv_breakpoints_read_support.tsv ${sample_id}_savana_run_breakpoints_read_support.tsv
        """
}

process runSavanaClassify {

    label 'savana'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(vcf_file), path(read_support_tsv_file), path(breakpoints_bedpe)
        path(custom_params_file)
        val(params_savana_classify)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_savana_classify_output.somatic.vcf"), path("${sample_id}_savana_classify_output.germline.vcf"), emit: f

    script:
        """
        savana classify \
            --vcf $vcf_file \
            --output ${sample_id}_savana_classify_output.vcf \
            --custom_params $custom_params_file \
            $params_savana_classify
        """
}
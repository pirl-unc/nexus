#!/usr/bin/env nextflow

process runExactoAnnotateVariants {

    label 'exacto_annotate_vars'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tsv_file)
        path(reference_gene_annotation_file)
        val(params_exacto_annotate_vars)
        val(output_dir)

    output:
        tuple val(sample_id), path("${tsv_file.baseName}_exacto_variant_annotation.tsv"), emit: f

    script:
        """
        exacto annotate-vars \
            --tsv-file $tsv_file \
            --reference-gene-annotation-file $reference_gene_annotation_file \
            --output-tsv-file ${tsv_file.baseName}_exacto_variant_annotation.tsv \
            --num-threads ${task.cpus} \
            $params_exacto_annotate_vars
        """
}

process runExactoCallSomaticDNAVariants {

    label 'exacto_call_dna_vars'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(tumor_sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        val(params_exacto_call_dna_vars)
        val(output_dir)

    output:
        tuple val(tumor_sample_id), path("${tumor_sample_id}_exacto_somatic_dna_variant_calls.tsv"), emit: f

    script:
        """
        exacto call-dna-vars \
            --bam-file $tumor_bam_file \
            --bam-bai-file $tumor_bam_bai_file \
            --control-bam-files $normal_bam_file \
            --control-bam-bai-files $normal_bam_bai_file \
            --mode case-specific \
            --output-tsv-file ${tumor_sample_id}_exacto_somatic_dna_variant_calls.tsv \
            --num-threads ${task.cpus} \
            $params_exacto_call_dna_vars
        """
}

process runExactoCallRNAVariants {

    label 'exacto_call_rna_vars'
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
        path(reference_gene_annotation_file)
        val(params_exacto_call_rna_vars)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${sample_id}_exacto_exons.tsv"),
              path("${sample_id}_exacto_read_filter_status.tsv"),
              path("${sample_id}_exacto_reference_transcript_matches.tsv"),
              path("${sample_id}_exacto_introns.tsv"),
              path("${sample_id}_exacto_transcripts_read_support.tsv"),
              path("${sample_id}_exacto_transcripts.tsv"),
              path("${sample_id}_exacto_rna_variant_calls.tsv"),
              emit: f

    script:
        """
        exacto call-rna-vars \
            --bam-file $bam_file \
            --bam-bai-file $bam_bai_file \
            --reference-genome-fasta-file $reference_genome_fasta_file \
            --reference-gene-annotation-file $reference_gene_annotation_file \
            --output-dir . \
            --output-prefix ${sample_id} \
            --num-threads ${task.cpus} \
            $params_exacto_call_rna_vars
        """
}

process runExactoIntegrateVariants {

    label 'exacto_integrate_vars'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), val(tumor_dna_id), val(tumor_rna_id), path(annotated_dna_variant_callset_tsv_file), path(rna_variant_callset_tsv_file)
        path(reference_gene_annotation_file)
        val(params_exacto_integrate_vars)
        val(output_dir)

    output:
        tuple val(sample_id), val(tumor_dna_id), val(tumor_rna_id), path("${tumor_dna_id}_${tumor_rna_id}_exacto_dna_rna_variant_integration.tsv"), emit: f

    script:
        """
        exacto integrate-vars \
            --annotated-dna-vars-tsv-file $annotated_dna_variant_callset_tsv_file \
            --rna-vars-tsv-file $rna_variant_callset_tsv_file \
            --reference-gene-annotation-file $reference_gene_annotation_file \
            --output-tsv-file ${tumor_dna_id}_${tumor_rna_id}_exacto_dna_rna_variant_integration.tsv \
            --num-threads ${task.cpus} \
            $params_exacto_integrate_vars
        """
}

process runExactoTranslateFastq {

    label 'exacto_translate'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_exacto_translate)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${sample_id}_exacto_translations.fasta.gz"),
              path("${sample_id}_exacto_translations.fasta.gz.fai"),
              path("${sample_id}_exacto_translations.fasta.gz.gzi"),
              path("${sample_id}_exacto_translations.tsv.gz"),
              emit: f

    script:
        """
        exacto translate \
            --fastq-file $fastq_file \
            --output-fasta-file ${sample_id}_exacto_translations.fasta.gz \
            --output-tsv-file ${sample_id}_exacto_translations.tsv.gz \
            --num-threads ${task.cpus} \
            --gzip yes \
            $params_exacto_translate
        """
}

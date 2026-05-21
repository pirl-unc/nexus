#!/usr/bin/env nextflow

process runExactoAnnotateVars {

    label 'exacto_annotate_vars'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/exacto/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tsv_file)
        path(reference_gene_annotation_file)
        val(reference_gene_annotation_source)
        val(reference_gene_annotation_assembly)
        val(reference_gene_annotation_version)
        val(params_exacto_annotate_vars)
        val(output_dir)

    output:
        tuple val(sample_id), path("${tsv_file.baseName}_exacto_variant_annotation.tsv"), emit: f

    script:
        """
        exacto annotate-vars \
            --tsv-file $tsv_file \
            --reference-gene-annotation-file $reference_gene_annotation_file \
            --reference-gene-annotation-source $reference_gene_annotation_source \
            --reference-gene-annotation-assembly $reference_gene_annotation_assembly \
            --reference-gene-annotation-version $reference_gene_annotation_version \
            --output-tsv-file ${tsv_file.baseName}_exacto_variant_annotation.tsv \
            --num-threads ${task.cpus} \
            $params_exacto_annotate_vars
        """
}

process runExactoBuildGenomeVarGraph {

    label 'exacto_build_genome_var_graph'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/exacto/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(variants_tsv_file)
        path(reference_genome_fasta_file)
        val(sequence_prefix)
        val(params_exacto_build_genome_var_graph)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_exacto_genome_var_graph.fasta"), emit: f

    script:
        """
        exacto build-genome-var-graph \
            --variants-tsv-file $variants_tsv_file \
            --fasta-file $reference_genome_fasta_file \
            --output-fasta-file ${sample_id}_exacto_genome_var_graph.fasta \
            --sequence-prefix $sequence_prefix \
            --num-threads ${task.cpus} \
            $params_exacto_build_genome_var_graph
        """
}

process runExactoBuildTranscriptomeVarGraph {

    label 'exacto_build_transcriptome_var_graph'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/exacto/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(transcript_structures_tsv_file)
        path(reference_genome_fasta_file)
        val(params_exacto_build_transcriptome_var_graph)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_exacto_transcriptome_var_graph.fasta"), emit: f

    script:
        """
        exacto build-transcriptome-var-graph \
            --transcript-structures-tsv-file $transcript_structures_tsv_file \
            --fasta-file $reference_genome_fasta_file \
            --output-fasta-file ${sample_id}_exacto_transcriptome_var_graph.fasta \
            --num-threads ${task.cpus} \
            $params_exacto_build_transcriptome_var_graph
        """
}

process runExactoCallGermlineDNAVars {

    label 'exacto_call_germline_dna_vars'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/exacto/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_gzi_file)
        val(params_exacto_call_germline_dna_vars)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_exacto_germline_dna_variant_calls.tsv"), emit: f

    script:
        """
        exacto call-germline-dna-vars \
            --bam-file $bam_file \
            --bam-bai-file $bam_bai_file \
            --fasta-file $reference_genome_fasta_file \
            --output-tsv-file ${sample_id}_exacto_germline_dna_variant_calls.tsv \
            --num-threads ${task.cpus} \
            $params_exacto_call_germline_dna_vars
        """
}

process runExactoCallSomaticDNAVars {

    label 'exacto_call_somatic_dna_vars'
    tag "${tumor_sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/exacto/",
        mode: 'copy'
    )

    input:
        tuple val(tumor_sample_id),
              path(tumor_bam_file), path(tumor_bam_bai_file),
              path(normal_bam_file), path(normal_bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_gzi_file)
        val(params_exacto_call_somatic_dna_vars)
        val(output_dir)

    output:
        tuple val(tumor_sample_id), path("${tumor_sample_id}_exacto_somatic_dna_variant_calls.tsv"), emit: f

    script:
        """
        exacto call-somatic-dna-vars \
            --bam-file $tumor_bam_file \
            --bam-bai-file $tumor_bam_bai_file \
            --control-bam-files $normal_bam_file \
            --control-bam-bai-files $normal_bam_bai_file \
            --fasta-file $reference_genome_fasta_file \
            --output-tsv-file ${tumor_sample_id}_exacto_somatic_dna_variant_calls.tsv \
            --num-threads ${task.cpus} \
            $params_exacto_call_somatic_dna_vars
        """
}

process runExactoCallPeptideVars {

    label 'exacto_call_peptide_vars'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/exacto/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(primary_structures_tsv_file)
        path(reference_proteome_fasta_file)
        val(params_exacto_call_peptide_vars)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${sample_id}_exacto_peptide_variant_calls.tsv"),
              path("${sample_id}_exacto_peptide_variant_calls.fasta"),
              emit: f

    script:
        """
        exacto call-peptide-vars \
            --primary-structures-tsv-file $primary_structures_tsv_file \
            --reference-fasta-file $reference_proteome_fasta_file \
            --output-tsv-file ${sample_id}_exacto_peptide_variant_calls.tsv \
            --output-fasta-file ${sample_id}_exacto_peptide_variant_calls.fasta \
            --num-threads ${task.cpus} \
            $params_exacto_call_peptide_vars
        """
}

process runExactoCallRNAVars {

    label 'exacto_call_rna_vars'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/exacto/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_gzi_file)
        path(reference_gene_annotation_file)
        val(reference_gene_annotation_source)
        val(reference_gene_annotation_assembly)
        val(reference_gene_annotation_version)
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
              path("${sample_id}_exacto_transcript_structures.tsv"),
              path("${sample_id}_exacto_rna_variant_calls.tsv"),
              emit: f

    script:
        """
        exacto call-rna-vars \
            --bam-file $bam_file \
            --bam-bai-file $bam_bai_file \
            --reference-genome-fasta-file $reference_genome_fasta_file \
            --reference-gene-annotation-file $reference_gene_annotation_file \
            --reference-gene-annotation-source $reference_gene_annotation_source \
            --reference-gene-annotation-assembly $reference_gene_annotation_assembly \
            --reference-gene-annotation-version $reference_gene_annotation_version \
            --output-dir . \
            --output-prefix ${sample_id} \
            --num-threads ${task.cpus} \
            $params_exacto_call_rna_vars
        """
}

process runExactoIntegrateVars {

    label 'exacto_integrate_vars'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/exacto/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id),
              val(tumor_dna_id), val(tumor_rna_id),
              path(annotated_dna_variant_callset_tsv_file),
              path(rna_variant_callset_tsv_file)
        path(reference_gene_annotation_file)
        val(reference_gene_annotation_source)
        val(reference_gene_annotation_assembly)
        val(reference_gene_annotation_version)
        val(params_exacto_integrate_vars)
        val(output_dir)

    output:
        tuple val(sample_id),
              val(tumor_dna_id), val(tumor_rna_id),
              path("${tumor_dna_id}_${tumor_rna_id}_exacto_dna_rna_variant_integration.tsv"),
              emit: f

    script:
        """
        exacto integrate-vars \
            --annotated-dna-vars-tsv-file $annotated_dna_variant_callset_tsv_file \
            --rna-vars-tsv-file $rna_variant_callset_tsv_file \
            --reference-gene-annotation-file $reference_gene_annotation_file \
            --reference-gene-annotation-source $reference_gene_annotation_source \
            --reference-gene-annotation-assembly $reference_gene_annotation_assembly \
            --reference-gene-annotation-version $reference_gene_annotation_version \
            --output-tsv-file ${tumor_dna_id}_${tumor_rna_id}_exacto_dna_rna_variant_integration.tsv \
            --num-threads ${task.cpus} \
            $params_exacto_integrate_vars
        """
}

process runExactoRemoveUnsplicedRNAs {

    label 'exacto_remove_unspliced_rnas'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/exacto/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(fasta_file)
        path(reference_gene_annotation_file)
        val(reference_gene_annotation_source)
        val(reference_gene_annotation_assembly)
        val(reference_gene_annotation_version)
        val(params_exacto_remove_unspliced_rnas)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${sample_id}_exacto_spliced.bam"),
              path("${sample_id}_exacto_spliced.bam.bai"),
              path("${sample_id}_exacto_spliced.fasta"),
              emit: f

    script:
        """
        exacto remove-unspliced-rnas \
            --bam-file $bam_file \
            --bam-bai-file $bam_bai_file \
            --fasta-file $fasta_file \
            --reference-gene-annotation-file $reference_gene_annotation_file \
            --reference-gene-annotation-source $reference_gene_annotation_source \
            --reference-gene-annotation-assembly $reference_gene_annotation_assembly \
            --reference-gene-annotation-version $reference_gene_annotation_version \
            --output-bam-file ${sample_id}_exacto_spliced.bam \
            --output-bam-bai-file ${sample_id}_exacto_spliced.bam.bai \
            --output-fasta-file ${sample_id}_exacto_spliced.fasta \
            --num-threads ${task.cpus} \
            $params_exacto_remove_unspliced_rnas
        """
}

process runExactoTranslateSeqs {

    label 'exacto_translate_seqs'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/exacto/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(input_file)
        val(input_type)                          // 'fastq' or 'fasta'
        val(strategy)                            // 'longest_orf' or 'all_orfs'
        val(params_exacto_translate_seqs)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${sample_id}_exacto_translations.fasta.gz"),
              path("${sample_id}_exacto_translations.tsv.gz"),
              emit: f

    script:
        def input_flag = (input_type == 'fastq') ? "--fastq-file" : "--fasta-file"
        """
        exacto translate-seqs \
            $input_flag $input_file \
            --strategy $strategy \
            --output-fasta-file ${sample_id}_exacto_translations.fasta.gz \
            --output-tsv-file ${sample_id}_exacto_translations.tsv.gz \
            --gzip yes \
            --num-threads ${task.cpus} \
            $params_exacto_translate_seqs
        """
}

process runExactoTranslateStructs {

    label 'exacto_translate_structs'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/exacto/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id),
              path(transcript_structures_tsv_file),
              path(rna_variant_calls_tsv_file),
              path(integrated_variants_tsv_file)
        val(strategy)                            // 'longest_orf' or 'all_orfs'
        val(params_exacto_translate_structs)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${sample_id}_exacto_primary_structures.tsv"),
              path("${sample_id}_exacto_primary_structures.fasta"),
              emit: f

    script:
        """
        exacto translate-structs \
            --transcript-structures-tsv-file $transcript_structures_tsv_file \
            --rna-variant-calls-tsv-file $rna_variant_calls_tsv_file \
            --integrated-variants-tsv-file $integrated_variants_tsv_file \
            --strategy $strategy \
            --output-tsv-file ${sample_id}_exacto_primary_structures.tsv \
            --output-fasta-file ${sample_id}_exacto_primary_structures.fasta \
            --num-threads ${task.cpus} \
            $params_exacto_translate_structs
        """
}

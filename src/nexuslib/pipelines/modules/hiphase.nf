#!/usr/bin/env nextflow

process runHiPhaseWith2VcfFiles {

    label 'hiphase'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id),
              path(bam_file),
              path(bam_bai_file),
              path(small_variants_vcf_gz_file),
              path(small_variants_vcf_gz_tbi_file),
              path(structural_variants_vcf_gz_file),
              path(structural_variants_vcf_gz_tbi_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${bam_file.baseName}_phased.bam"),
              path("${bam_file.baseName}_phased.bam.bai"),
              path("${small_variants_vcf_gz_file.baseName}_phased.vcf.gz"),
              path("${structural_variants_vcf_gz_file.baseName}_phased.vcf.gz"),
              path("${sample_id}_hiphase_summary.tsv"),
              path("${sample_id}_hiphase_stats.tsv"),
              path("${sample_id}_hiphase_blocks.tsv"),
              path("${sample_id}_hiphase_haplotag.tsv"), emit: f

    script:
        """
        hiphase \
            --bam $bam_file \
            --output-bam ${bam_file.baseName}_phased.bam \
            --vcf $small_variants_vcf_gz_file \
            --output-vcf ${small_variants_vcf_gz_file.baseName}_phased.vcf.gz \
            --vcf $structural_variants_vcf_gz_file \
            --output-vcf ${structural_variants_vcf_gz_file.baseName}_phased.vcf.gz \
            --reference $reference_genome_fasta_file \
            --sample-name $sample_id \
            --summary-file ${sample_id}_hiphase_summary.tsv \
            --stats-file ${sample_id}_hiphase_stats.tsv \
            --blocks-file ${sample_id}_hiphase_blocks.tsv \
            --haplotag-file ${sample_id}_hiphase_haplotag.tsv \
            --phase-singletons \
            --threads ${task.cpus}
        """
}

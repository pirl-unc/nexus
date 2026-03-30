#!/usr/bin/env nextflow

process runHiPhaseWith2VcfFiles {

    label 'hiphase'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
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
        tuple val(sample_id), path("${sample_id}_hiphase_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_hiphase_outputs/
        hiphase \
            --bam $bam_file \
            --output-bam ${sample_id}_hiphase_outputs/${bam_file.baseName}_phased.bam \
            --vcf $small_variants_vcf_gz_file \
            --output-vcf ${sample_id}_hiphase_outputs/${small_variants_vcf_gz_file.baseName}_phased.vcf.gz \
            --vcf $structural_variants_vcf_gz_file \
            --output-vcf ${sample_id}_hiphase_outputs/${structural_variants_vcf_gz_file.baseName}_phased.vcf.gz \
            --reference $reference_genome_fasta_file \
            --sample-name $sample_id \
            --summary-file ${sample_id}_hiphase_outputs/${sample_id}_hiphase_summary.tsv \
            --stats-file ${sample_id}_hiphase_outputs/${sample_id}_hiphase_stats.tsv \
            --blocks-file ${sample_id}_hiphase_outputs/${sample_id}_hiphase_blocks.tsv \
            --haplotag-file ${sample_id}_hiphase_outputs/${sample_id}_hiphase_haplotag.tsv \
            --phase-singletons \
            --threads ${task.cpus}
        """
}

#!/usr/bin/env nextflow

process runHiPhaseWith2VcfFiles {

    label 'hiphase'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "*.vcf.gz"
    )
    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "*.vcf.gz.tbi"
    )

    // HiPhase TSV files (haplotag, summary, stats, blocks) — always published
    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "*.tsv"
    )

    // Haplotagged BAM + index — published only when mode includes 'bam'
    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "*_phased.bam",
        enabled: ((params.haplotag_output ?: 'bam').toString().toLowerCase() in ['bam', 'both'])
    )
    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "*_phased.bam.bai",
        enabled: ((params.haplotag_output ?: 'bam').toString().toLowerCase() in ['bam', 'both'])
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
        val(subdir)

    output:
        tuple val(sample_id),
              path("${small_variants_vcf_gz_file.baseName}_phased.vcf.gz"),
              path("${small_variants_vcf_gz_file.baseName}_phased.vcf.gz.tbi"),
              path("${structural_variants_vcf_gz_file.baseName}_phased.vcf.gz"),
              path("${structural_variants_vcf_gz_file.baseName}_phased.vcf.gz.tbi"),
              path("${bam_file.baseName}_phased.bam"),
              path("${bam_file.baseName}_phased.bam.bai"),
              path("${sample_id}_hiphase_summary.tsv"),
              path("${sample_id}_hiphase_stats.tsv"),
              path("${sample_id}_hiphase_blocks.tsv"),
              path("${sample_id}_hiphase_haplotag.tsv"), emit: f

    script:
        """
        # HiPhase requires --sample-name to match the VCF sample column.
        # DeepVariant/pbsv copy the SM tag from the BAM's read group, which
        # may differ from the workflow's sample_id (e.g. SM:hg007 vs sample_id
        # hg007_dna_pacbio-nist_minimap2). Extract the actual name from the
        # small-variants VCF header so HiPhase can find it.
        VCF_SAMPLE=\$(zcat $small_variants_vcf_gz_file | grep -m1 '^#CHROM' | awk '{print \$10}')
        echo "Using VCF sample name: \$VCF_SAMPLE"

        hiphase \
            --bam $bam_file \
            --output-bam ${bam_file.baseName}_phased.bam \
            --vcf $small_variants_vcf_gz_file \
            --output-vcf ${small_variants_vcf_gz_file.baseName}_phased.vcf.gz \
            --vcf $structural_variants_vcf_gz_file \
            --output-vcf ${structural_variants_vcf_gz_file.baseName}_phased.vcf.gz \
            --reference $reference_genome_fasta_file \
            --sample-name \$VCF_SAMPLE \
            --summary-file ${sample_id}_hiphase_summary.tsv \
            --stats-file ${sample_id}_hiphase_stats.tsv \
            --blocks-file ${sample_id}_hiphase_blocks.tsv \
            --haplotag-file ${sample_id}_hiphase_haplotag.tsv \
            --phase-singletons \
            --threads ${task.cpus}
        """
}

#!/usr/bin/env nextflow

/*
 * HapCUT2 phaser. Runs in two stages that are tightly coupled:
 *   1. `extractHAIRS` — extracts haplotype-relevant read fragments from the BAM
 *   2. `HAPCUT2`     — phases the variants using those fragments
 * The fragment file is intermediate; the published outputs are the phased
 * haplotype block file and the phased VCF.
 */

process runHapCUT2 {

    label 'hapcut2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${sample_id}_hapcut2_haplotypes.txt"
    )

    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${sample_id}_hapcut2_haplotypes.txt.phased.vcf.gz"
    )

    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${sample_id}_hapcut2_haplotypes.txt.phased.vcf.gz.tbi"
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(vcf_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(read_technology)             // 'pacbio' | 'ont' | 'illumina'
        val(params_extracthairs)
        val(params_hapcut2)
        val(output_dir)
        val(subdir)

    output:
        tuple val(sample_id),
              path("${sample_id}_hapcut2_haplotypes.txt"),
              path("${sample_id}_hapcut2_haplotypes.txt.phased.vcf.gz"),
              path("${sample_id}_hapcut2_haplotypes.txt.phased.vcf.gz.tbi"),
              emit: f

    script:
        // Pick the right read-technology flag for extractHAIRS.
        def tech = read_technology?.toString()?.toLowerCase() ?: 'pacbio'
        def tech_flag = ''
        if (tech == 'pacbio')        tech_flag = '--pacbio 1'
        else if (tech == 'ont')      tech_flag = '--ont 1'
        else if (tech == 'illumina') tech_flag = ''                // illumina is the default; no flag
        else error "ERROR: read_technology must be one of [pacbio, ont, illumina] (got: '${read_technology}')."

        """
        # ----------------------------------------------------------
        # 0. Normalize the input VCF to plain .vcf
        #    (extractHAIRS / HAPCUT2 expect uncompressed VCF).
        # ----------------------------------------------------------
        if [[ "${vcf_file}" == *.gz ]]; then
            zcat ${vcf_file} > ${sample_id}_input.vcf
        else
            cp ${vcf_file} ${sample_id}_input.vcf
        fi

        # ----------------------------------------------------------
        # 1. extractHAIRS — produce per-read fragment file
        # ----------------------------------------------------------
        extractHAIRS \
            ${tech_flag} \
            --bam ${bam_file} \
            --VCF ${sample_id}_input.vcf \
            --ref ${reference_genome_fasta_file} \
            --out ${sample_id}_hapcut2_fragments.txt \
            ${params_extracthairs}

        # ----------------------------------------------------------
        # 2. HAPCUT2 — phase variants using the fragment file.
        #    --outvcf 1 produces a *.phased.VCF alongside the haplotype
        #    block file (named after the --output value).
        # ----------------------------------------------------------
        HAPCUT2 \
            --fragments ${sample_id}_hapcut2_fragments.txt \
            --VCF ${sample_id}_input.vcf \
            --output ${sample_id}_hapcut2_haplotypes.txt \
            --outvcf 1 \
            ${params_hapcut2}

        # ----------------------------------------------------------
        # 3. bgzip + tabix-index the phased VCF.
        #    HAPCUT2 writes uppercase .phased.VCF; bgzip-c straight to the
        #    lowercase .phased.vcf.gz target avoids the case-insensitive
        #    filesystem rename trap on macOS bind mounts.
        # ----------------------------------------------------------
        bgzip -c ${sample_id}_hapcut2_haplotypes.txt.phased.VCF \
            > ${sample_id}_hapcut2_haplotypes.txt.phased.vcf.gz
        tabix -f -p vcf ${sample_id}_hapcut2_haplotypes.txt.phased.vcf.gz
        """
}

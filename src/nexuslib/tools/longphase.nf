#!/usr/bin/env nextflow

process runLongphaseHaplotag {

    label 'longphase'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${bam_file.baseName}_haplotagged.bam",
        enabled: ((params.haplotag_output ?: 'bam').toString().toLowerCase() in ['bam', 'both'])
    )

    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${bam_file.baseName}_haplotagged.bam.bai",
        enabled: ((params.haplotag_output ?: 'bam').toString().toLowerCase() in ['bam', 'both'])
    )

    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${bam_file.baseName}_haplotagged_haplotag.tsv.gz",
        enabled: ((params.haplotag_output ?: 'bam').toString().toLowerCase() in ['tsv', 'both'])
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(phased_vcf_file), path(phased_vcf_tbi_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_longphase)
        val(output_dir)
        val(subdir)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_haplotagged.bam"), path("${bam_file.baseName}_haplotagged.bam.bai"), emit: f
        path("${bam_file.baseName}_haplotagged_haplotag.tsv.gz"), emit: tsv

    script:
        // LongPhase has no native haplotag-list output flag, so we extract a
        // 2-column (readname, HP) TSV from the haplotagged BAM via samtools.
        // `longphase haplotag -o <prefix>` writes <prefix>.bam.
        """
        longphase haplotag \
            --reference $reference_genome_fasta_file \
            -s $phased_vcf_file \
            -b $bam_file \
            -t ${task.cpus} \
            -o ${bam_file.baseName}_haplotagged \
            $params_longphase
        samtools index -@ ${task.cpus} -b ${bam_file.baseName}_haplotagged.bam ${bam_file.baseName}_haplotagged.bam.bai

        # Extract a (readname, HP) TSV from the haplotagged BAM.
        # 2 columns; '.' for reads without an HP tag.
        {
            printf "# readname\\thaplotype\\n"
            samtools view ${bam_file.baseName}_haplotagged.bam | awk -v OFS='\\t' '{
                hp="."
                for (i=12; i<=NF; i++) if (\$i ~ /^HP:i:/) { hp=substr(\$i, 6); break }
                print \$1, hp
            }'
        } | gzip > ${bam_file.baseName}_haplotagged_haplotag.tsv.gz
        """
}

process runLongphasePhase {

    label 'longphase'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${sample_id}_longphase_phased.vcf.gz"
    )

    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${sample_id}_longphase_phased.vcf.gz.tbi"
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(vcf_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_longphase)
        val(output_dir)
        val(subdir)

    output:
        tuple val(sample_id), path("${sample_id}_longphase_phased.vcf.gz"), path("${sample_id}_longphase_phased.vcf.gz.tbi"), emit: f

    script:
        // LongPhase reads plain VCFs; some builds don't accept .gz directly,
        // so decompress to a sibling file when needed.
        """
        if [[ "${vcf_file}" == *.gz ]]; then
            gunzip -c ${vcf_file} > input.vcf
            in_vcf=input.vcf
        else
            in_vcf=${vcf_file}
        fi

        longphase phase \
            --reference $reference_genome_fasta_file \
            -s \${in_vcf} \
            -b $bam_file \
            -t ${task.cpus} \
            -o ${sample_id}_longphase_phased \
            $params_longphase

        # `longphase phase -o <prefix>` writes <prefix>.vcf (uncompressed).
        bgzip ${sample_id}_longphase_phased.vcf
        tabix -p vcf ${sample_id}_longphase_phased.vcf.gz
        """
}

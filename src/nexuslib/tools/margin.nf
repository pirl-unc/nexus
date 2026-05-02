#!/usr/bin/env nextflow

process runMarginPhase {

    label 'margin'
    tag "${sample_id}"
    debug true

    // Phased VCF + index
    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${sample_id}_margin_phased.vcf.gz"
    )
    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${sample_id}_margin_phased.vcf.gz.tbi"
    )

    // Haplotagged BAM + index (gated by --haplotag_output, like other tools).
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
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(vcf_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(margin_params_json_file)
        val(params_margin)
        val(output_dir)
        val(subdir)

    output:
        tuple val(sample_id),
              path("${sample_id}_margin_phased.vcf.gz"),
              path("${sample_id}_margin_phased.vcf.gz.tbi"),                                emit: phased_vcf
        tuple val(sample_id),
              path("${bam_file.baseName}_haplotagged.bam"),
              path("${bam_file.baseName}_haplotagged.bam.bai"),                             emit: haplotagged_bam
        path("${bam_file.baseName}_haplotagged_haplotag.tsv.gz"),                           emit: tsv

    script:
        // Margin v2.3.1 dropped the standalone `haplotag` subcommand:
        //   margin phase <BAM> <REF.fa> <VCF> <params.json> -o <prefix>
        // By default, `phase` emits BOTH the phased VCF (<prefix>.phased.vcf)
        // and the haplotagged BAM (<prefix>.haplotagged.bam). Use
        // `-V --skipPhasedVCF` or `-M --skipHaplotypeBAM` to suppress either.
        """
        # Margin needs the input VCF uncompressed.
        if [[ "${vcf_file}" == *.gz ]]; then
            gunzip -c ${vcf_file} > input.vcf
            in_vcf=input.vcf
        else
            in_vcf=${vcf_file}
        fi

        # Margin errors out when the input VCF has zero entries. That is a
        # legitimate situation (e.g. a region with no variants, or a caller
        # that produced an empty VCF on this sample). Detect it early and emit
        # minimal valid outputs so the rest of the workflow can proceed.
        n_entries=\$(grep -cv '^#' \${in_vcf} || true)
        if [ "\${n_entries:-0}" -lt 1 ]; then
            echo "WARNING: input VCF \${in_vcf} has 0 entries — emitting empty phased VCF and untagged BAM copy."
            # Empty phased VCF (header-only) + index.
            grep '^#' \${in_vcf} > ${sample_id}_margin_phased.vcf || true
            if [ ! -s ${sample_id}_margin_phased.vcf ]; then
                printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n' \
                    > ${sample_id}_margin_phased.vcf
            fi
            bgzip ${sample_id}_margin_phased.vcf
            tabix -p vcf ${sample_id}_margin_phased.vcf.gz

            # Haplotagged BAM: copy input BAM (no HP tags). Sort + index.
            samtools sort -@ ${task.cpus} \
                -o ${bam_file.baseName}_haplotagged.bam \
                ${bam_file}
            samtools index -@ ${task.cpus} -b \
                ${bam_file.baseName}_haplotagged.bam \
                ${bam_file.baseName}_haplotagged.bam.bai

            # Empty haplotag TSV (header only).
            printf "# readname\\thaplotype\\n" | gzip \
                > ${bam_file.baseName}_haplotagged_haplotag.tsv.gz
            exit 0
        fi

        margin phase \
            ${bam_file} \
            ${reference_genome_fasta_file} \
            \${in_vcf} \
            ${margin_params_json_file} \
            -t ${task.cpus} \
            -o ${sample_id}_margin \
            ${params_margin}

        # Phased VCF: `margin phase -o <prefix>` writes <prefix>.phased.vcf.
        mv ${sample_id}_margin.phased.vcf ${sample_id}_margin_phased.vcf
        bgzip ${sample_id}_margin_phased.vcf
        tabix -p vcf ${sample_id}_margin_phased.vcf.gz

        # Haplotagged BAM: margin v2.3.1's `phase` writes <prefix>.haplotagged.bam
        # by default. Rename to align with the `<bam_basename>_haplotagged.bam`
        # convention used by whatshap/longphase so downstream consumers can find
        # it by a stable name. Margin emits an unsorted BAM; sort before indexing.
        samtools sort -@ ${task.cpus} \
            -o ${bam_file.baseName}_haplotagged.bam \
            ${sample_id}_margin.haplotagged.bam
        samtools index -@ ${task.cpus} -b \
            ${bam_file.baseName}_haplotagged.bam \
            ${bam_file.baseName}_haplotagged.bam.bai

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

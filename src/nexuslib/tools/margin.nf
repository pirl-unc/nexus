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
            printf "# readname\\thaplotype\\tphaseset\\n" | gzip \
                > ${bam_file.baseName}_haplotagged_haplotag.tsv.gz
            exit 0
        fi

        # Margin requires each read ID to have AT MOST one primary alignment.
        # If two primaries share a read name, margin polishes through to the
        # merge step and then fails with:
        #   "Expected three tokens in header line, got 2
        #    This usually means you have multiple primary alignments with the
        #    same read ID."
        # ...wasting hours of polishing. Detect duplicates up front and, if
        # present, write a deduplicated BAM (keep the first primary per read;
        # preserve all secondary / supplementary / unmapped records) and feed
        # that to margin. Other haplotaggers (whatshap, longphase) are not
        # sensitive to this so we keep the fix local to runMarginPhase rather
        # than touching the aligner.
        margin_bam=${bam_file}
        n_dup=\$(samtools view -F 0x904 ${bam_file} | cut -f1 | sort | uniq -d | wc -l)
        if [ "\${n_dup}" -gt 0 ]; then
            echo "WARNING: \${n_dup} read IDs have multiple primary alignments. Deduplicating BAM before margin."
            {
                samtools view -H ${bam_file}
                # Primary mapped: keep first occurrence per read ID.
                samtools view -F 0x904 ${bam_file} | awk '!seen[\$1]++'
                # Unmapped primaries (have 0x4, no 0x100/0x800).
                samtools view -f 0x004 -F 0x900 ${bam_file}
                # Secondary alignments.
                samtools view -f 0x100 ${bam_file}
                # Supplementary alignments.
                samtools view -f 0x800 ${bam_file}
            } | samtools view -bS - | samtools sort -@ ${task.cpus} -o margin_input.bam -
            samtools index -@ ${task.cpus} margin_input.bam
            margin_bam=margin_input.bam
        fi

        margin phase \
            \${margin_bam} \
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

        # Extract a (readname, HP, PS) TSV from the haplotagged BAM.
        # 3 columns; '.' for reads without an HP or PS tag.
        {
            printf "# readname\\thaplotype\\tphaseset\\n"
            samtools view ${bam_file.baseName}_haplotagged.bam | awk -v OFS='\\t' '{
                hp="."; ps="."
                for (i=12; i<=NF; i++) {
                    if      (\$i ~ /^HP:i:/) hp=substr(\$i, 6)
                    else if (\$i ~ /^PS:i:/) ps=substr(\$i, 6)
                }
                print \$1, hp, ps
            }'
        } | gzip > ${bam_file.baseName}_haplotagged_haplotag.tsv.gz
        """
}

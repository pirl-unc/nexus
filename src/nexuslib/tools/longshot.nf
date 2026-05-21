#!/usr/bin/env nextflow

process runLongshot {

    label 'longshot'
    tag "${sample_id}"
    debug true

    publishDir(  // VCF — always
        path: "${output_dir}/${sample_id}/longshot/",
        mode: 'copy',
        pattern: "${sample_id}_longshot.vcf"
    )
    publishDir(  // BAM — only when mode includes 'bam'
        path: "${output_dir}/${sample_id}/longshot/",
        mode: 'copy',
        pattern: "${sample_id}_longshot.bam",
        enabled: ((params.haplotag_output ?: 'bam').toString().toLowerCase() in ['bam', 'both'])
    )
    publishDir(  // BAI — only when mode includes 'bam'
        path: "${output_dir}/${sample_id}/longshot/",
        mode: 'copy',
        pattern: "${sample_id}_longshot.bam.bai",
        enabled: ((params.haplotag_output ?: 'bam').toString().toLowerCase() in ['bam', 'both'])
    )
    publishDir(  // haplotag TSV — only when mode includes 'tsv'
        path: "${output_dir}/${sample_id}/longshot/",
        mode: 'copy',
        pattern: "${sample_id}_longshot_haplotag.tsv.gz",
        enabled: ((params.haplotag_output ?: 'bam').toString().toLowerCase() in ['tsv', 'both'])
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_longshot)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_longshot.vcf"), path("${sample_id}_longshot.bam"), path("${sample_id}_longshot.bam.bai"), emit: f
        path("${sample_id}_longshot_haplotag.tsv.gz"), optional: true, emit: tsv

    script:
        // Only generate the TSV when at least one TSV-publishing mode is active.
        def emit_tsv = (params.haplotag_output ?: 'bam').toString().toLowerCase() in ['tsv', 'both']
        """
        longshot \
            --bam $bam_file \
            --ref $reference_genome_fasta_file \
            --sample_id $sample_id \
            --out ${sample_id}_longshot.vcf \
            --out_bam ${sample_id}_longshot.bam \
            $params_longshot

        # If longshot finds zero variants it skips writing the output BAM,
        # which breaks the output: declaration. Fall back to the input BAM
        # in that case so downstream consumers always have a BAM to point at.
        if [ ! -f ${sample_id}_longshot.bam ]; then
            cp $bam_file ${sample_id}_longshot.bam
        fi
        samtools index -@ ${task.cpus} -b ${sample_id}_longshot.bam ${sample_id}_longshot.bam.bai

        # Extract a (readname, HP, PS) TSV from the phased BAM when requested.
        # 3 columns; '.' for reads without an HP or PS tag.
        if [ "${emit_tsv}" = "true" ]; then
            {
                printf "# readname\\thaplotype\\tphaseset\\n"
                samtools view ${sample_id}_longshot.bam | awk -v OFS='\\t' '{
                    hp="."; ps="."
                    for (i=12; i<=NF; i++) {
                        if      (\$i ~ /^HP:i:/) hp=substr(\$i, 6)
                        else if (\$i ~ /^PS:i:/) ps=substr(\$i, 6)
                    }
                    print \$1, hp, ps
                }'
            } | gzip > ${sample_id}_longshot_haplotag.tsv.gz
        fi
        """
}

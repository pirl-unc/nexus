#!/usr/bin/env nextflow

process runLongcallR {

    label 'longcallr'
    tag "${sample_id}"
    debug true

    publishDir(  // VCF — always
        path: "${output_dir}/${sample_id}/",
        mode: 'copy',
        pattern: "${sample_id}_longcallr.vcf"
    )
    publishDir(  // BAM — only when mode includes 'bam'
        path: "${output_dir}/${sample_id}/",
        mode: 'copy',
        pattern: "${sample_id}_longcallr.phased.bam",
        enabled: ((params.haplotag_output ?: 'bam').toString().toLowerCase() in ['bam', 'both'])
    )
    publishDir(  // haplotag TSV — only when mode includes 'tsv'
        path: "${output_dir}/${sample_id}/",
        mode: 'copy',
        pattern: "${sample_id}_longcallr_haplotag.tsv.gz",
        enabled: ((params.haplotag_output ?: 'bam').toString().toLowerCase() in ['tsv', 'both'])
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genes_gtf_file)
        val(preset)
        val(params_longcallr)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_longcallr.phased.bam"), path("${sample_id}_longcallr.vcf"), emit: f
        path("${sample_id}_longcallr_haplotag.tsv.gz"), optional: true, emit: tsv

    script:
        def emit_tsv = (params.haplotag_output ?: 'bam').toString().toLowerCase() in ['tsv', 'both']
        """
        longcallR \
            --bam-path $bam_file \
            --ref-path $reference_genome_fasta_file \
            --annotation $reference_genes_gtf_file \
            --preset $preset \
            --output ${sample_id}_longcallr \
            --threads ${task.cpus} \
            $params_longcallr

        # Extract a (readname, HP, PS) TSV from the phased BAM when requested.
        # 3 columns; '.' for reads without an HP or PS tag.
        if [ "${emit_tsv}" = "true" ]; then
            {
                printf "# readname\\thaplotype\\tphaseset\\n"
                samtools view ${sample_id}_longcallr.phased.bam | awk -v OFS='\\t' '{
                    hp="."; ps="."
                    for (i=12; i<=NF; i++) {
                        if      (\$i ~ /^HP:i:/) hp=substr(\$i, 6)
                        else if (\$i ~ /^PS:i:/) ps=substr(\$i, 6)
                    }
                    print \$1, hp, ps
                }'
            } | gzip > ${sample_id}_longcallr_haplotag.tsv.gz
        fi
        """
}

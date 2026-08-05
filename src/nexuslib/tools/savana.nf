#!/usr/bin/env nextflow

process runSavanaRun {

    label 'savana'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(contigs_file)
        val(params_savana_run)
        val(output_dir)

    output:
        tuple val(sample_id), path("savana/"), emit: f

    script:
        """
        mkdir -p savana/
        savana run \
            -t $tumor_bam_file \
            -n $normal_bam_file \
            --ref $reference_genome_fasta_file \
            --ref_index $reference_genome_fasta_fai_file \
            --contigs $contigs_file \
            --threads ${task.cpus} \
            --outdir savana/ \
            --sample $sample_id \
            $params_savana_run
        """
}

process runSavanaClassify {

    label 'savana'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/savana/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(savana_run_dir)
        path(custom_params_file)
        val(params_savana_classify)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_savana_classify_output.somatic.vcf"), path("${sample_id}_savana_classify_output.germline.vcf"), emit: f

    script:
        """
        breakpoints_vcf=${savana_run_dir}/${sample_id}.sv_breakpoints.vcf

        n_breakpoints=0
        if [ -f "\${breakpoints_vcf}" ]; then
            n_breakpoints=\$(grep -cv '^#' "\${breakpoints_vcf}" || true)
        else
            echo "WARNING: \${breakpoints_vcf} not found — 'savana run' emitted no breakpoints VCF for ${sample_id}."
        fi

        if [ "\${n_breakpoints:-0}" -lt 1 ]; then
            echo "WARNING: 'savana run' produced 0 SV breakpoints for ${sample_id} — skipping 'savana classify' and emitting empty somatic/germline VCFs."
            for out in ${sample_id}_savana_classify_output.somatic.vcf \
                       ${sample_id}_savana_classify_output.germline.vcf; do
                if [ -f "\${breakpoints_vcf}" ]; then
                    # Reuse the breakpoints header so contig/INFO/FORMAT definitions
                    # survive, tagging why the file is empty.
                    awk '/^#CHROM/ { print "##nexus_note=savana_classify_skipped_zero_breakpoints" } { print }' \
                        "\${breakpoints_vcf}" > "\${out}"
                fi
                if [ ! -s "\${out}" ]; then
                    printf '##fileformat=VCFv4.2\\n##nexus_note=savana_classify_skipped_zero_breakpoints\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n' \
                        > "\${out}"
                fi
            done
            exit 0
        fi

        savana classify \
            --vcf "\${breakpoints_vcf}" \
            --output ${sample_id}_savana_classify_output.vcf \
            --custom_params $custom_params_file \
            $params_savana_classify
        """
}
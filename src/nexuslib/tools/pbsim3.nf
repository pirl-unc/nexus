#!/usr/bin/env nextflow

process runPbsim3DNA {

    label 'pbsim3_dna'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fasta_file)
        path(pbsim3_model_file)
        val(params_pbsim3_mode)
        val(params_pbsim3)
        val(min_contig_length)
        val(output_dir)

    output:
        // All per-contig BAMs
        tuple val(sample_id),
              path("*_pbsim3*.bam"),
              emit: bam

        // All per-contig MAF.GZs
        tuple val(sample_id),
              path("*_pbsim3*.maf.gz"),
              emit: maf_gz

        // All per-contig REF files
        tuple val(sample_id),
              path("*_pbsim3*.ref"),
              emit: ref

        // Final merged BAM
        tuple val(sample_id),
              path("${sample_id}_pbsim3.merged.bam"),
              emit: merged_bam

        // Final merged FASTQ
        tuple val(sample_id),
              path("${sample_id}_pbsim3.merged.fq.gz"),
              emit: merged_fastq

    script:
        """
        mkdir -p contigs/

        seqkit seq -m ${min_contig_length} ${fasta_file} | seqkit split -i -O contigs

        ls -1 contigs/*.fasta | \
          xargs -n 1 -P ${task.cpus} -I {} bash -lc '
            fa="{}"
            contig=\$(basename "\$fa" .fa)

            pbsim \
              --method ${params_pbsim3_mode} \
              --${params_pbsim3_mode} ${pbsim3_model_file} \
              --genome "\$fa" \
              --id-prefix "${sample_id}_\${contig}_" \
              --prefix "${sample_id}_\${contig}_pbsim3" \
              ${params_pbsim3}
          '

        if [[ "${params_pbsim3_mode}" == "errhmm" ]]; then
            samtools merge -@ ${task.cpus} -o ${sample_id}_pbsim3.merged.bam ${sample_id}*pbsim3*.bam
            : > ${sample_id}_pbsim3.merged.fq.gz
        else
            cat ${sample_id}*.fq.gz > ${sample_id}_pbsim3.merged.fq.gz
            : > ${sample_id}_pbsim3.merged.bam
        fi
        """
}

process runPbsim3RNA {

    label 'pbsim3_rna'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(transcript_file)
        path(pbsim3_model_file)
        val(params_pbsim3_mode)
        val(params_pbsim3)
        val(output_dir)

    output:
        // Per-chunk BAMs (errhmm only)
        tuple val(sample_id),
              path("${sample_id}_chunk_*_pbsim3*.bam"),
              emit: bam,             optional: true

        // Per-chunk MAFs (errhmm only)
        tuple val(sample_id),
              path("${sample_id}_chunk_*_pbsim3*.maf.gz"),
              emit: maf_gz,          optional: true

        // Final merged BAM
        tuple val(sample_id),
              path("${sample_id}_pbsim3.merged.bam"),
              emit: merged_bam

        // Final merged FASTQ (qshmm only; empty for errhmm)
        tuple val(sample_id),
              path("${sample_id}_pbsim3.merged.fq.gz"),
              emit: merged_fastq

    script:
        """
        mkdir -p chunks/

        # Split the transcript file (one transcript per line) evenly across CPUs
        total=\$(wc -l < ${transcript_file})
        per=\$(( (total + ${task.cpus} - 1) / ${task.cpus} ))
        split -d -a 4 -l \$per ${transcript_file} chunks/chunk_

        ls -1 chunks/chunk_* | \
          xargs -n 1 -P ${task.cpus} -I {} bash -lc '
            chunk="{}"
            name=\$(basename "\$chunk")

            pbsim \
              --method ${params_pbsim3_mode} \
              --${params_pbsim3_mode} ${pbsim3_model_file} \
              --transcript "\$chunk" \
              --id-prefix "${sample_id}_\${name}" \
              --prefix "${sample_id}_\${name}_pbsim3" \
              ${params_pbsim3}
          '

        # The ajslee/pbsim3 container emits BAM+maf.gz when pass-num>1 (multi-pass
        # CCS) and fq.gz+maf.gz when pass-num=1. Detect at runtime so both work.
        if compgen -G "${sample_id}_chunk_*_pbsim3*.bam" > /dev/null; then
            samtools merge -@ ${task.cpus} -o ${sample_id}_pbsim3.merged.bam ${sample_id}_chunk_*_pbsim3*.bam
            : > ${sample_id}_pbsim3.merged.fq.gz
        else
            # per-chunk files are already gzipped; concatenated gzip is a valid gzip stream
            cat ${sample_id}_chunk_*_pbsim3*.fq.gz > ${sample_id}_pbsim3.merged.fq.gz
            : > ${sample_id}_pbsim3.merged.bam
        fi
        """
}

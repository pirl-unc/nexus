#!/usr/bin/env nextflow

process runPbsim3DNA {

    label 'pbsim3_dna'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}",
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
        path: "${output_dir}/${sample_id}",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(transcript_file)
        path(pbsim3_model_file)
        val(params_pbsim3_mode)
        val(params_pbsim3)
        val(output_dir)

   output:
        tuple val(sample_id), path("${sample_id}_pbsim3.bam"), path("${sample_id}_pbsim3.maf.gz"), emit: f

   script:
        """
        pbsim \
            --method ${params_pbsim3_mode} \
            --${params_pbsim3_mode} ${pbsim3_model_file} \
            --transcript $transcript_file \
            --id-prefix ${sample_id} \
            --prefix "${sample_id}_pbsim3" \
            ${params_pbsim3}
        """
}

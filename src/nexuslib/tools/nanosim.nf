#!/usr/bin/env nextflow

process runNanoSimGenome {

    label 'nanosim'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/nanosim/",
        mode: 'copy'
    )

    input:
        // Stage the input FASTA under input/ so that the renamed FASTA
        // we write at the work-dir root (${sample_id}_reference_genome.fa)
        // can never collide with the input filename (which would otherwise
        // truncate the user's source file via the staged symlink).
        tuple val(sample_id), path(fasta_file, stageAs: 'input/*')
        val(model_prefix)
        val(params_nanosim)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_nanosim_genome_aligned_reads.fastq.gz"), path("${sample_id}_nanosim_genome_aligned_error_profile"), emit: f

    script:
        """
        RENAMED=_nanosim_ref_genome_renamed_${sample_id}.fa
        awk -v sid="${sample_id}" '
            /^>/ {
                sub(/^>/, "")
                split(\$0, a, /[ \t]/)
                print ">" a[1] "_" sid
                next
            }
            { print }
        ' $fasta_file > "\$RENAMED"

        simulator.py genome \
            --ref_g "\$RENAMED" \
            --output ${sample_id}_nanosim_genome \
            --model_prefix $model_prefix \
            --fastq \
            --num_threads ${task.cpus} \
            $params_nanosim

        gzip ${sample_id}_nanosim_genome_aligned_reads.fastq
        """
}

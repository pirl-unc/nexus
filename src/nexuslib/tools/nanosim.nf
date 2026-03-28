#!/usr/bin/env nextflow

process runNanoSimGenome {

    label 'nanosim'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fasta_file)
        val(model_prefix)
        val(params_nanosim)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_nanosim_genome_aligned_reads.fastq.gz"), path("${sample_id}_nanosim_genome_aligned_reads_renamed.fastq.gz"), path("${sample_id}_nanosim_genome_read_names_map.tsv"), path("${sample_id}_nanosim_genome_aligned_error_profile"), emit: f

    script:
        """
        simulator.py genome \
            --ref_g $fasta_file \
            --output ${sample_id}_nanosim_genome \
            --model_prefix $model_prefix \
            --fastq \
            --num_threads ${task.cpus} \
            $params_nanosim

        gzip ${sample_id}_nanosim_genome_aligned_reads.fastq

        # Rename read IDs
        PREFIX=${sample_id}
        IN=${sample_id}_nanosim_genome_aligned_reads.fastq.gz
        OUT=${sample_id}_nanosim_genome_aligned_reads_renamed.fastq.gz
        MAP="${sample_id}_nanosim_genome_read_names_map.tsv"

        printf "original_read_name\tnew_read_name\n" >> "\$MAP"

        zcat "\$IN" | awk -v prefix="\$PREFIX" -v map="\$MAP" '
        BEGIN { i=1 }
        NR%4==1 {
            orig=\$0
            sub(/^@/, "", orig)
            split(orig, a, /[ \t]/)
            old=a[1]
            new=prefix "_" i
            print "@" new
            print old "\t" new >> map
            i++
            next
        }
        { print }
        ' | gzip -c > "\$OUT"
        """
}

#!/usr/bin/env nextflow

process runNanoCaller {

    label 'nanocaller'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_nanocaller)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_nanocaller_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_nanocaller_outputs/
        NanoCaller \
            --bam $bam_file \
            --ref $reference_genome_fasta_file \
            --output ${sample_id}_nanocaller_outputs/ \
            --cpu ${task.cpus} \
            --prefix ${sample_id}_nanocaller \
            --sample ${sample_id} \
            $params_nanocaller
        """
}

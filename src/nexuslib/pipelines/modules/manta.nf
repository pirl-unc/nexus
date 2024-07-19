#!/usr/bin/env nextflow

process runMantaTumorNormal {

    label 'manta'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_manta_config)
        val(params_manta_run)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_manta_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_manta_outputs/
        configManta.py \
            --tumorBam $tumor_bam_file \
            --normalBam $normal_bam_file \
            --referenceFasta $reference_genome_fasta_file \
            --runDir ${sample_id}_manta_outputs/ \
            $params_manta_config
        ${sample_id}_manta_outputs/runWorkflow.py \
            -j ${task.cpus} \
            -g ${task.memory.toGiga()} \
            $params_manta_run
        """
}

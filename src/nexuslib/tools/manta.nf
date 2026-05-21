#!/usr/bin/env nextflow

process runMantaSomatic {

    label 'manta'
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
        val(params_manta_config)
        val(params_manta_run)
        val(output_dir)

    output:
        tuple val(sample_id), path("manta/"), emit: f

    script:
        """
        mkdir -p manta/
        configManta.py \
            --tumorBam $tumor_bam_file \
            --normalBam $normal_bam_file \
            --referenceFasta $reference_genome_fasta_file \
            --runDir manta/ \
            $params_manta_config
        manta/runWorkflow.py \
            -j ${task.cpus} \
            -g ${task.memory.toGiga()} \
            $params_manta_run
        """
}

process runMantaGermline {

    label 'manta'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_manta_config)
        val(params_manta_run)
        val(output_dir)

    output:
        tuple val(sample_id), path("manta/"), emit: f

    script:
        """
        mkdir -p manta/
        configManta.py \
            --bam $bam_file \
            --referenceFasta $reference_genome_fasta_file \
            --runDir manta/ \
            $params_manta_config
        manta/runWorkflow.py \
            -j ${task.cpus} \
            -g ${task.memory.toGiga()} \
            $params_manta_run
        """
}

#!/usr/bin/env nextflow

process runStrelka2GermlineMode {

    label 'strelka2'
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
        path(reference_genome_fasta_gzi_file)
        val(params_strelka2)
        val(output_dir)

   output:
         tuple val(sample_id), path("strelka2/"), emit: f

   script:
        """
        mkdir -p strelka2/
        configureStrelkaGermlineWorkflow.py \
            --bam ${bam_file} \
            --referenceFasta ${reference_genome_fasta_file} \
            --runDir strelka2/ \
            $params_strelka2
        python strelka2/runWorkflow.py -m local -j ${task.cpus}
        """
}

process runStrelka2SomaticMode {

    label 'strelka2'
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
        path(reference_genome_fasta_gzi_file)
        val(params_strelka2)
        val(output_dir)

    output:
         tuple val(sample_id), path("strelka2/"), emit: f

    script:
        """
        mkdir -p strelka2/
        configureStrelkaSomaticWorkflow.py \
            --normalBam $normal_bam_file \
            --tumorBam $tumor_bam_file \
            --referenceFasta $reference_genome_fasta_file \
            --runDir strelka2/ \
            $params_strelka2
        python strelka2/runWorkflow.py -m local -j ${task.cpus}
        """
}

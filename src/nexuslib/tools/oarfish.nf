#!/usr/bin/env nextflow

process runOarfishFastqMode {

    label 'oarfish'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        path(reference_transcriptome_fasta_file)
        val(params_oarfish)
        val(output_dir)

    output:
        tuple val(sample_id), path("oarfish/"), emit: f

    script:
        """
        mkdir -p oarfish/
        oarfish \
            --reference $reference_transcriptome_fasta_file \
            --threads ${task.cpus} \
            $params_oarfish \
            --output oarfish/${sample_id}_oarfish \
            --reads $fastq_file
        """
}

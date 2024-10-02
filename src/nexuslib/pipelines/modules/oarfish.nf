#!/usr/bin/env nextflow

process runOarfishRawReadMode {

    label 'oarfish'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        path(reference_transcriptome_fasta_file)
        val(params_oarfish)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_oarfish_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_oarfish_outputs/
        oarfish \
            --reference $reference_transcriptome_fasta_file \
            --threads ${task.cpus} \
            $params_oarfish \
            --output ${sample_id}_oarfish_outputs/${sample_id}_oarfish \
            --reads $fastq_file
        """
}

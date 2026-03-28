#!/usr/bin/env nextflow

process runPindel {

    label 'pindel'
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
        val(params_pindel)
        val(output_dir)

    output:
        tuple val(sample_id), emit: f

    script:
        """
        echo "${bam_file} 500 ${sample_id}" > pindel.input
        pindel \
            -f $reference_genome_fasta_file \
            -i pindel.input \
            -o ${sample_id}_pindel \
            -T ${task.cpus} \
            $params_pindel
        """
}

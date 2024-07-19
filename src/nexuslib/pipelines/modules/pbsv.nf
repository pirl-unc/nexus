#!/usr/bin/env nextflow

process runPbsv {

    label 'pbsv'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(reference_genome_fasta_file)
        val(params_pbsv_discover)
        val(params_pbsv_call)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_pbsv.vcf"), emit: f

    script:
        """
        pbsv discover \
            $params_pbsv_discover \
            $bam_file \
            ${sample_id}_pbsv.svsig.gz
        pbsv call \
            --num-threads ${task.cpus} \
            $params_pbsv_call \
            $reference_genome_fasta_file \
            ${sample_id}_pbsv.svsig.gz \
            ${sample_id}_pbsv.vcf
        """
}
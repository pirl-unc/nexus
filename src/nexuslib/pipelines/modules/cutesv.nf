#!/usr/bin/env nextflow

process runCuteSV {

    label 'cutesv'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(reference_genome_fasta_file)
        val(cutesv)
        val(cutesv_params)
        val(output_dir)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_cutesv.vcf"), emit: f

    script:
        def cutesv_params_ = cutesv_params == true ? '' : cutesv_params

        """
        mkdir -p output/

        $cutesv \
            --threads ${task.cpus} \
            --sample $sample_id \
            $cutesv_params_ \
            $bam_file \
            $reference_genome_fasta_file \
            ${bam_file.baseName}_cutesv.vcf \
            output/
        """
}
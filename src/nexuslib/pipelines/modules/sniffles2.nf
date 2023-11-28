#!/usr/bin/env nextflow

process runSniffles2 {

    label 'sniffles2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(reference_genome_fasta_file)
        val(sniffles2)
        val(sniffles2_params)
        val(output_dir)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_sniffles2.vcf"), emit: f

    script:
        def sniffles2_params_ = sniffles2_params == true ? '' : sniffles2_params

        """
        $sniffles2 \
            --input $bam_file \
            --vcf ${bam_file.baseName}_sniffles2.vcf \
            --snf ${bam_file.baseName}_sniffles2.snf \
            --sample-id ${sample_id} \
            --reference ${reference_genome_fasta_file} \
            --threads ${task.cpus} \
            $sniffles2_params_
        """
}

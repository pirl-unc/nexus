#!/usr/bin/env nextflow

process runClair3RNA {

    label 'clair3rna'
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
        val(params_clair3rna)
        val(output_dir)

    output:
        tuple val(sample_id), path("clair3rna/"), emit: f

    script:
        """
        mkdir -p clair3rna/
        run_clair3_rna \
            --bam_fn $bam_file \
            --ref_fn $reference_genome_fasta_file \
            --output_dir clair3rna/ \
            --threads ${task.cpus} \
            $params_clair3rna
        """
}

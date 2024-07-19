#!/usr/bin/env nextflow

process runArriba {

    label 'arriba'
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
        path(gtf_file)
        path(protein_domains_gff3_file)
        val(params_arriba)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_arriba.tsv"), emit: f

    script:
        """
        arriba \
            -x $bam_file \
            -g $gtf_file \
            -a $reference_genome_fasta_file \
            -p $protein_domains_gff3_file \
            -o ${sample_id}_arriba.tsv \
            $params_arriba
        """
}
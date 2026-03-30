#!/usr/bin/env nextflow

process runIsotools {

    label 'isotools'
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
        path(reference_genes_gtf_file)
        val(params_isotools)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_isotools_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_isotools_outputs/
        printf "sample_name\tfile_name\tgroup\n" >> sample.tsv
        printf "%s\t%s\t%s\n" "$sample_id" "$bam_file" "$sample_id" >> sample.tsv

        run_isotools \
            --anno $reference_genes_gtf_file \
            --genome $reference_genome_fasta_file \
            --samples sample.tsv \
            --file_prefix ${sample_id}_isotools_outputs/${sample_id} \
            --gtf_out
            $params_isotools
        """
}

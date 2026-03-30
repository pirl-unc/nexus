#!/usr/bin/env nextflow

process runTalon {

    label 'talon'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genes_gtf_file)
        val(params_talon_initdb)
        val(params_talon)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_talon_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_talon_outputs/

        talon_initialize_database \
            --f $reference_genes_gtf_file \
            --a reference_annotation \
            --g reference_genome \
            --o ${sample_id}_talon_outputs/reference \
            $params_talon_initdb

        printf "%s,%s,%s,%s\n" "$sample_id" "$sample_id" "$sample_id" "$bam_file" >> config.csv

        talon \
            --f config.csv \
            --db ${sample_id}_talon_outputs/reference.db \
            --build reference_genome \
            --threads ${task.cpus} \
            --o ${sample_id}_talon_outputs/${sample_id} \
            $params_talon
        """
}

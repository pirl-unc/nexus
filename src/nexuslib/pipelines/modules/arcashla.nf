#!/usr/bin/env nextflow

process runArcasHlaPairedEndMode {

    label 'arcashla'
    tag "${sample_id}"
    debug true

   publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), val(bam_bai_file)
        val(output_dir)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_arcashla_genotype.json"), emit: f

    script:

        """
        mkdir -p output/
        arcasHLA extract $bam_file -o output/ -t ${task.cpus} -v
        arcasHLA genotype output/*.extracted.1.fq.gz output/*.extracted.2.fq.gz --genes all -o output/ -t ${task.cpus} -v
        cp output/*.genotype.json ${bam_file.baseName}_arcashla_genotype.json
        """
}
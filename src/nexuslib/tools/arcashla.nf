#!/usr/bin/env nextflow

process runArcasHlaPairedEndMode {

    label 'arcashla'
    tag "${sample_id}"
    debug true

   publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), val(bam_bai_file)
        val(output_dir)

    output:
        tuple val(sample_id), path("arclashla/"), emit: f

    script:

        """
        mkdir -p arclashla/
        arcasHLA extract $bam_file -o arclashla/ -t ${task.cpus} -v
        arcasHLA genotype arclashla/*.extracted.1.fq.gz arclashla/*.extracted.2.fq.gz --genes all -o arclashla/ -t ${task.cpus} -v
        """
}
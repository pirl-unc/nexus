#!/usr/bin/env nextflow

process runStringTie3 {

    label 'stringtie3'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(gtf_file)
        val(params_stringtie3)
        val(output_dir)

    output:
        tuple val(sample_id), emit: f

    script:
        """
        stringtie \
            -o ${sample_id}_stringtie3.gtf \
            -G $gtf_file \
            -p ${task.cpus} \
            $params_stringtie3 \
            $bam_file
        """
}

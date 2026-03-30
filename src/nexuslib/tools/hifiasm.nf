#!/usr/bin/env nextflow

process runHifiasm {

    label 'hifiasm'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_files)
        val(params_hifiasm)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_hifiasm_outputs/"), emit: f

    script:
        """
        READS="${fastq_files.join(' ')}"
        mkdir -p ${sample_id}_hifiasm_outputs/
        hifiasm \
            $params_hifiasm \
            -t ${task.cpus} \
            -o ${sample_id}_hifiasm_outputs/${sample_id}_hifiasm.asm \
            \$READS
        """
}

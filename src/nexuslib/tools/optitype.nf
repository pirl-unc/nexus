#!/usr/bin/env nextflow

process runOptiType {

    label 'optitype'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(params_optitype)
        val(output_dir)

    output:
        tuple val(sample_id), path("optitype/"), emit: f

    script:
        """
        mkdir -p optitype/

        sed '/^\\[mapping\\]/,/^\\[/ s/^threads=.*/threads=${task.cpus}/' \
            /usr/local/bin/OptiType/config.ini > config.ini

        OptiTypePipeline.py \
            --input $fastq_file_1 \
            --input $fastq_file_2 \
            --outdir optitype/ \
            --config config.ini \
            $params_optitype
        """
}

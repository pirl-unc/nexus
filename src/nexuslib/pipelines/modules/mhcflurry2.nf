#!/usr/bin/env nextflow

process runMhcflurry2 {

    label 'mhcflurry2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(mhcflurry2_input_csv_file)
        val(params_mhcflurry2_predict)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_mhcflurry2_outputs.csv"), emit: f

    script:
        """
        mhcflurry-predict \
            --out ${sample_id}_mhcflurry2_outputs.csv \
            --models /opt/mhcflurry2/models/ \
            $params_mhcflurry2_predict \
            $mhcflurry2_input_csv_file
        """
}

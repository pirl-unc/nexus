#!/usr/bin/env nextflow

process runKallistoQuantSingleEndReads {

    label 'kallisto'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        path(kallisto_index_file)
        val(params_kallisto_quant_fragment_length)
        val(params_kallisto_quant_sd)
        val(params_kallisto_quant)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_kallisto_output/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_kallisto_output/
        kallisto quant \
            --index=$kallisto_index_file \
            --output-dir=${sample_id}_kallisto_output/ \
            --fragment-length=$params_kallisto_quant_fragment_length \
            --sd=$params_kallisto_quant_sd \
            --single \
            --threads=${task.cpus} \
            $params_kallisto_quant \
            $fastq_file
        """
}

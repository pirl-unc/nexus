#!/usr/bin/env nextflow

process runIsoquant {

    label 'isoquant'
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        val(reference_genome_fasta_file)
        val(isoquant)
        val(isoquant_params)
        val(output_dir)

    output:
        path('isoquant_outputs/'), emit: f

    script:
        def isoquant_params_ = isoquant_params == true ? '' : isoquant_params

        """
        mkdir -p isoquant_outputs/
        $isoquant \
          --reference $reference_genome_fasta_file \
          --threads ${task.cpus} \
          -o isoquant_outputs/ \
          $isoquant_params_
        """
}


#!/usr/bin/env nextflow

process runTranSigner {

    label 'transigner'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        path(reference_transcriptome_fasta_file)
        val(params_transigner_align)
        val(params_transigner_prefilter)
        val(params_transigner_em)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_transigner_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_transigner_outputs/

        transigner align \
            -q $fastq_file \
            -t $reference_transcriptome_fasta_file \
            -d ${sample_id}_transigner_outputs/ \
            -o ${sample_id}_alignment.bam \
            -p ${task.cpus} \
            $params_transigner_align

        transigner prefilter \
            -a ${sample_id}_alignment.bam \
            -t $reference_transcriptome_fasta_file \
            -o ${sample_id}_transigner_outputs/ \
            $params_transigner_prefilter

        transigner em \
            -s ${sample_id}_transigner_outputs/scores.tsv \
            -i ${sample_id}_transigner_outputs/ti.pkl \
            -o ${sample_id}_transigner_outputs/ \
            $params_transigner_em
        """
}

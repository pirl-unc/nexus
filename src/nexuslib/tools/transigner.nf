#!/usr/bin/env nextflow

process runTranSigner {

    label 'transigner'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
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
        tuple val(sample_id), path("transigner/"), emit: f

    script:
        """
        mkdir -p transigner/

        transigner align \
            -q $fastq_file \
            -t $reference_transcriptome_fasta_file \
            -d transigner/ \
            -o ${sample_id}_alignment.bam \
            -p ${task.cpus} \
            $params_transigner_align

        transigner prefilter \
            -a ${sample_id}_alignment.bam \
            -t $reference_transcriptome_fasta_file \
            -o transigner/ \
            $params_transigner_prefilter

        transigner em \
            -s transigner/scores.tsv \
            -i transigner/ti.pkl \
            -o transigner/ \
            $params_transigner_em
        """
}

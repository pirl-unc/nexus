#!/usr/bin/env nextflow

process runIsonClust3 {

    label 'isonclust3'
    tag "${sample_id}"
    debug true
    // Stage inputs as real copies, not symlinks.
    stageInMode 'copy'

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_isonclust3)
        val(output_dir)

    output:
        tuple val(sample_id), path("isonclust3/"), emit: f

    script:
        """
        mkdir -p isonclust3/
        gunzip -c $fastq_file > ${sample_id}.fastq
        isONclust3 \
            --fastq \$PWD/${sample_id}.fastq \
            --outfolder \${PWD}/isonclust3/ \
            $params_isonclust3
        """
}

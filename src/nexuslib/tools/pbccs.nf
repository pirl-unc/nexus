#!/usr/bin/env nextflow

process runPbccs {

    label 'pbccs'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/pbccs/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(subreads_bam_file)
        val(params_pbccs)
        val(output_dir)

    output:
        tuple val(sample_id), path("${subreads_bam_file.baseName}.ccs.bam"), path("${subreads_bam_file.baseName}.ccs.fastq.gz"), emit: f

    script:
        """
        pbindex --num-threads ${task.cpus} $subreads_bam_file

        ccs \
            --num-threads ${task.cpus} \
            $params_pbccs \
            $subreads_bam_file \
            ${subreads_bam_file.baseName}.ccs.bam

        samtools fastq --threads \$(( ${task.cpus} / 2 )) ${subreads_bam_file.baseName}.ccs.bam \
            | pigz -p \$(( ${task.cpus} / 2 )) > ${subreads_bam_file.baseName}.ccs.fastq.gz
        """
}
#!/usr/bin/env nextflow

process runPbccs {

    label 'pbccs'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
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

        for i in {1..20}; do
            ccs \
                --num-threads ${task.cpus} \
                --chunk \${i}/20 \
                $params_pbccs \
                $subreads_bam_file \
                ${subreads_bam_file.baseName}.ccs.\${i}.bam
        done

        samtools cat ${subreads_bam_file.baseName}.ccs.*.bam \
            | tee ${subreads_bam_file.baseName}.ccs.bam \
            | samtools fastq --threads \$(( ${task.cpus} / 2 )) - \
            | pigz -p \$(( ${task.cpus} / 2 )) > ${subreads_bam_file.baseName}.ccs.fastq.gz
        """
}
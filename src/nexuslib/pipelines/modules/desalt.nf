#!/usr/bin/env nextflow

process runDesalt {

    label 'desalt'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file)
        val(desalt)
        val(desalt_index)
        val(desalt_info_file)
        val(desalt_params)
        val(samtools)

    output:
        tuple val(sample_id), path("${sample_id}_desalt_sorted.bam"), path("${sample_id}_desalt_sorted.bam.bai"), emit: f

    script:
        """
        $desalt aln \
            --thread ${task.cpus} \
            -G $desalt_info_file \
            $desalt_params \
            -f ${sample_id}_desalt_temp \
            --output ${sample_id}_desalt.sam \
            $desalt_index \
            $fastq_file
        $samtools sort -@ {task.samtools_cpus} -m ${task.samtools_memory.toGiga()}G -O bam -o ${sample_id}_desalt_sorted.bam ${sample_id}_desalt.sam
        $samtools index -b ${sample_id}_desalt_sorted.bam ${sample_id}_desalt_sorted.bam.bai
        """
}

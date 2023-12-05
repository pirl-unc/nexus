#!/usr/bin/env nextflow

process runSamtoolsSamToBam {

    label 'samtools_sort'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(sam_file)
        val(samtools)

    output:
        tuple val(sample_id), path("${sam_file.baseName}.bam"), path("${sam_file.baseName}.bam.bai"), emit: f

    script:
        """
        $samtools sort -@ ${task.cpus} -m ${task.samtools_memory.toGiga()}G -O bam -o ${sam_file.baseName}.bam $sam_file
        $samtools index -b ${sam_file.baseName}.bam ${sam_file.baseName}.bam.bai
        """
}

process runSamtoolsSort {

    label 'samtools_sort'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(samtools)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_sorted.bam"), path("${bam_file.baseName}_sorted.bam.bai"), emit: f

    script:
        """
        $samtools sort -@ ${task.cpus} -m ${task.samtools_memory.toGiga()}G -O bam -o ${bam_file.baseName}_sorted.bam $bam_file
        $samtools index -b ${bam_file.baseName}_sorted.bam ${bam_file.baseName}_sorted.bam.bai
        """
}

process runSamtoolsIndex {

    label 'samtools_index'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file)
        val(samtools)

    output:
        tuple val(sample_id), path("${bam_file}"), path("${bam_file.baseName}.bam.bai"), emit: f

    script:
        """
        $samtools index -b ${bam_file} ${bam_file.baseName}.bam.bai
        """
}

process runSamtoolsCalmd {

    label 'samtools_calmd'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(samtools)
        val(reference_genome_fasta_file)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_mdtagged.bam"), path("${bam_file.baseName}_mdtagged.bam.bai"), emit: f

    script:
        """
        $samtools calmd -b --threads ${task.cpus} $bam_file $reference_genome_fasta_file > ${bam_file.baseName}_mdtagged.bam
        $samtools index -b ${bam_file.baseName}_mdtagged.bam ${bam_file.baseName}_mdtagged.bam.bai
        """
}

process runSamtoolsMarkdup {

    label 'samtools_markdup'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file)
        val(samtools)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_markeddup.bam"), path("${bam_file.baseName}_markeddup.bam.bai"), emit: f

    script:
        """
        $samtools sort -@ ${task.cpus} -m ${task.samtools_memory.toGiga()}G $bam_file -o ${bam_file.baseName}_position_sorted.bam
        $samtools markdup --threads ${task.cpus} ${bam_file.baseName}_position_sorted.bam ${bam_file.baseName}_markeddup.bam
        $samtools index -b ${bam_file.baseName}_markeddup.bam ${bam_file.baseName}_markeddup.bam.bai
        """
}

process runSamtoolsFixmate {

    label 'samtools_fixmate'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(samtools)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_fixmate.bam"), emit: f

    script:
        """
        $samtools sort -n -@ ${task.cpus} -m ${task.samtools_memory.toGiga()}G $bam_file -o ${bam_file.baseName}_name_sorted.bam
        $samtools fixmate -m --threads ${task.cpus} ${bam_file.baseName}_name_sorted.bam ${bam_file.baseName}_fixmate.bam
        """
}

process runSamtoolsCoverage {

    label 'samtools_coverage'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(samtools)
        val(min_mapping_quality)
        val(min_base_quality)
        val(output_dir)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_samtools_coverage.txt"), emit: f

    script:
        """
        $samtools coverage \
            --min-MQ ${min_mapping_quality} \
            --min-BQ ${min_base_quality} \
            $bam_file \
            --output ${bam_file.baseName}_samtools_coverage.txt
        """
}

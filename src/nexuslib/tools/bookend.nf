#!/usr/bin/env nextflow

process runBookendLabelSingleEnd {

    label 'bookend'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/bookend/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_bookend_label)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_bookend_label.fq"), emit: f

    script:
        """
        bookend label \
            --single_out ${sample_id}_bookend_label.fq \
            $params_bookend_label \
            $fastq_file
        """
}

process runBookendELR {

    label 'bookend'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/bookend/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_bookend_elr)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_bookend.elr"), emit: f

    script:
        """
        bookend elr \
            --output ${sample_id}_bookend \
            --source $sample_id \
            --genome $reference_genome_fasta_file \
            $params_bookend_elr \
            $bam_file
        """
}

process runBookendCondense {

    label 'bookend'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/bookend/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(elr_file)
        val(params_bookend_condense)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_bookend_condense.elr"), emit: f

    script:
        """
        bookend condense \
            --output ${sample_id}_bookend_condense.elr \
            $params_bookend_condense \
            $elr_file
        """
}

process runBookendAssemble {

    label 'bookend'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/bookend/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(params_bookend_assemble)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_bookend_assembly.gtf"), emit: f

    script:
        """
        bookend assemble \
            --source ${sample_id} \
            --output ${sample_id}_bookend_assembly.gtf \
            $params_bookend_assemble \
            $bam_file
        """
}

process runBookendFasta {

    label 'bookend'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/bookend/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(gtf_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_bookend_fasta)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_bookend_assembly.fasta"), emit: f

    script:
        """
        bookend fasta \
            --output ${sample_id}_bookend_assembly.fasta \
            --genome $reference_genome_fasta_file \
            $params_bookend_fasta \
            $gtf_file
        """
}

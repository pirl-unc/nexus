#!/usr/bin/env nextflow

process runBwaMem2 {

    label 'bwamem2'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_dict_file)
        path(reference_genome_fasta_0123_file)
        path(reference_genome_fasta_amb_file)
        path(reference_genome_fasta_ann_file)
        path(reference_genome_fasta_bwt_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_pac_file)
        val(platform_tag)
        val(platform_unit_tag)
        val(library_tag)

    output:
        tuple val(sample_id), path("${sample_id}.sam"), emit: f

    script:
        """
        bwa-mem2 mem -t ${task.cpus} \
            -R "@RG\tID:$sample_id\tSM:$sample_id\tPL:$platform_tag\tLB:$library_tag\tPU:$platform_unit_tag" \
            $reference_genome_fasta_file \
            $fastq_file_1 \
            $fastq_file_2 > ${sample_id}.sam
        """
}

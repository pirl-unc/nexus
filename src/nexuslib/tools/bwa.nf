#!/usr/bin/env nextflow

process runBwaIndex {

    label 'bwa_index'
    debug true

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.amb"), emit: amb_file
        path("${fasta_file}.ann"), emit: ann_file
        path("${fasta_file}.bwt"), emit: bwt_file
        path("${fasta_file}.pac"), emit: pac_file
        path("${fasta_file}.sa"), emit: sa_file

    script:
        """
        bwa index $fasta_file
        """
}
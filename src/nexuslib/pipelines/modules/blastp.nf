#!/usr/bin/env nextflow

process runBlastp {

    label 'blast'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fasta_file)
        path(reference_proteome_fasta_file)
        val(params_blastp)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_blastp_output.txt"), emit: f

    script:
        """
        if file "$reference_proteome_fasta_file" | grep -q 'gzip'; then
            gunzip -c ${reference_proteome_fasta_file} > reference.fasta
        else
            cp ${reference_proteome_fasta_file} reference.fasta
        fi
        makeblastdb -in reference.fasta -dbtype prot -out blastpdb
        blastp \
            -query $fasta_file \
            -db blastpdb \
            -out ${sample_id}_blastp_output.txt \
            -outfmt 7 \
            -num_threads ${task.cpus} \
            $params_blastp
        """
}

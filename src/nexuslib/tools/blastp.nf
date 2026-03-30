#!/usr/bin/env nextflow

process runBlastp {

    label 'blast'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fasta_file)
        path(blastdb_dir)
        val(blastdb_name)
        val(params_blastp)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_blastp_output.txt"), emit: f

    script:
        def blastdb_dir_str = blastdb_dir.toString()
        def blastdb_basename = blastdb_dir_str.tokenize(File.separator).last()
        def is_gzipped = fasta_file.name.endsWith(".gz")
        def query_command = is_gzipped ? "zcat $fasta_file" : "cat $fasta_file"

        """
        $query_command | blastp \
            -query - \
            -db $blastdb_basename/$blastdb_name \
            -out ${sample_id}_blastp_output.txt \
            -outfmt 7 \
            -num_threads ${task.cpus} \
            $params_blastp
        """
}

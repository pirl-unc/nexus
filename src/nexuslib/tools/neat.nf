#!/usr/bin/env nextflow

process runNeat {

    label 'neat'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fasta_file)
        val(params_neat)
        val(min_contig_length)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_neat_read*.fq.gz"), emit: f

    script:
        """
        mkdir -p contigs/

        seqkit seq -m ${min_contig_length} ${fasta_file} | seqkit split -i -O contigs

        ls -1 contigs/*.fasta | \
          xargs -n 1 -P ${task.cpus} -I {} bash -lc '
            fa="{}"
            contig=\$(basename "\$fa" .fa)

            /opt/conda/bin/python /opt/neat/gen_reads.py \
                -r "\$fa" \
                -o "${sample_id}_\${contig}" \
                ${params_neat}
          '

        cat ${sample_id}_*read1.fq.gz > ${sample_id}_neat_read1.fq.gz
        cat ${sample_id}_*read2.fq.gz > ${sample_id}_neat_read2.fq.gz
        """
}

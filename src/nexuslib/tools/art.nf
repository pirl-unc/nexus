#!/usr/bin/env nextflow

process runArtIlluminaPairedEnd {

    label 'art'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fasta_file)
        val(params_art_illumina)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_art-illumina.merged.r1.fq.gz"), path("${sample_id}_art-illumina.merged.r2.fq.gz"), emit: f

    script:
        """
        mkdir -p contigs/

        seqkit seq ${fasta_file} | seqkit split -i -O contigs

        ls -1 contigs/*.fasta | \
          xargs -n 1 -P ${task.cpus} -I {} bash -lc '
            fa="{}"
            contig=\$(basename "\$fa" .fa)

            art_illumina \
              $params_art_illumina \
              -p \
              -i "\$fa" \
              -d ${sample_id}- \
              -o ${sample_id}_\${contig}_art-illumina.r

            gzip ${sample_id}_\${contig}_art-illumina.r1.fq
            gzip ${sample_id}_\${contig}_art-illumina.r2.fq
          '

        cat ${sample_id}*.r1.fq.gz > ${sample_id}_art-illumina.merged.r1.fq.gz
        cat ${sample_id}*.r2.fq.gz > ${sample_id}_art-illumina.merged.r2.fq.gz
        """
}

#!/usr/bin/env nextflow

process runGeluster {

    label 'geluster'
    tag "${sample_id}"
    debug true
    // Stage inputs as real copies, not symlinks.
    stageInMode 'copy'

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_geluster)
        val(output_dir)

    output:
        tuple val(sample_id), path("geluster/"), emit: f
        tuple val(sample_id), path("geluster/${sample_id}_geluster_clusters.tsv"), emit: tsv

    script:
        """
        gunzip -c $fastq_file > ${sample_id}_long_read.fastq
        GeLuster \
            --reads ${sample_id}_long_read.fastq \
            --threads ${task.cpus} \
            --output_dir \${PWD}/geluster \
            $params_geluster

        # GeLuster.tsv lines are '<fastq header> ,gene_cluster_<N>', and reads that
        # cluster alone are written to GeLuster_singleton.tsv with no cluster
        # number. Merge both into the 'cluster_id', 'read_name' TSV that the
        # assembly_*_clustered subworkflows expect, giving every singleton its own
        # cluster id after the numbered ones. read_name is the first whitespace
        # token without '@'/'>', which is what pysam reports as FastxFile
        # entry.name.
        awk 'BEGIN { OFS="\\t"; print "cluster_id", "read_name"; max_id = -1 }
             { name = \$1; sub(/^[@>]/, "", name) }
             FILENAME ~ /singleton/ && NF { print ++max_id, name; next }
             match(\$0, /,gene_cluster_[0-9]+\$/) {
                 id = substr(\$0, RSTART + 14) + 0
                 if (id > max_id) max_id = id
                 print id, name
             }' geluster/GeLuster.tsv geluster/GeLuster_singleton.tsv \
            > geluster/${sample_id}_geluster_clusters.tsv
        """
}

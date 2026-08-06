#!/usr/bin/env nextflow

process runIsonClust3 {

    label 'isonclust3'
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
        val(params_isonclust3)
        val(output_dir)

    output:
        tuple val(sample_id), path("isonclust3/"), emit: f
        tuple val(sample_id), path("isonclust3/${sample_id}_isonclust3_clusters.tsv"), emit: tsv

    script:
        """
        mkdir -p isonclust3/
        gunzip -c $fastq_file > ${sample_id}.fastq
        isONclust3 \
            --fastq \$PWD/${sample_id}.fastq \
            --outfolder \${PWD}/isonclust3/ \
            $params_isonclust3

        # isONclust3 writes '<outfolder>/clustering/final_clusters.tsv' as
        # '<cluster_id>\\t<read_id>' with no header. Add the 'cluster_id',
        # 'read_name' header the assembly_*_clustered subworkflows expect.
        awk 'BEGIN { OFS="\\t"; print "cluster_id", "read_name" } NF >= 2 { print \$1, \$2 }' \
            isonclust3/clustering/final_clusters.tsv \
            > isonclust3/${sample_id}_isonclust3_clusters.tsv
        """
}

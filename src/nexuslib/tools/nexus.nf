#!/usr/bin/env nextflow

process runNexusFilterRNABloom2Transcripts {

    label 'nexus_filter_rnabloom2_transcripts'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/nexus/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(rnabloom2_output_dir)
        val(params_nexus_filter_rnabloom2_transcripts)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${sample_id}_rnabloom2_filtered_reads.tsv"),
              path("${sample_id}_rnabloom2_filtered_transcripts.tsv"),
              path("${sample_id}_rnabloom2_filtered.fasta"),
              path("${sample_id}_rnabloom2_filtered.fastq.gz"),
              emit: f

    script:
        """
        nexus_filter_rnabloom2_transcripts \
            --assembly4-pol-fasta-file ${rnabloom2_output_dir}/rnabloom.longreads.assembly4.pol.fa \
            --assembly3-map-paf-file ${rnabloom2_output_dir}/rnabloom.longreads.assembly3.map.paf.gz \
            --output-reads-tsv-file ${sample_id}_rnabloom2_filtered_reads.tsv \
            --output-transcripts-tsv-file ${sample_id}_rnabloom2_filtered_transcripts.tsv \
            --output-fasta-file ${sample_id}_rnabloom2_filtered.fasta \
            --output-fastq-file ${sample_id}_rnabloom2_filtered.fastq.gz \
            $params_nexus_filter_rnabloom2_transcripts
        """
}
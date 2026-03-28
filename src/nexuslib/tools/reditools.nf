#!/usr/bin/env nextflow

process runReditoolsDenovo {

    label 'reditools_denovo'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_gtf_file)
        val(params_reditoolsdenovo)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_reditools_outputs/"), emit: f

    script:
        """
        REDItoolDenovo.py \
            $params_reditoolsdenovo \
            -i $bam_file \
            -f $reference_genome_fasta_file \
            -o reditools_outputs \
            -X $reference_gtf_file \
            -t ${task.cpus}
        original_dir=\$(find reditools_outputs/ -maxdepth 1 -type d -name 'denovo_*' | head -n 1)
        mv "\$original_dir" "${sample_id}_reditools_outputs"
        """
}

process runReditoolsAnnotateTable {

    label 'reditools_annotate'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tsv_file)
        path(reference_gtf_file)
        val(params_reditools_annotatetable)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_reditools_annotated.tsv"), emit: f

    script:
        """
        AnnotateTable.py \
            -i $tsv_file \
            -a $reference_gtf_file \
            -o ${sample_id}_reditools_annotated.tsv \
            $params_reditools_annotatetable
        """
}

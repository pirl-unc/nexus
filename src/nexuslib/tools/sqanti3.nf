#!/usr/bin/env nextflow

process runSqanti3FastaMode {

    label 'sqanti3'
    tag "${sample_id}_${method}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fasta_file)
        path(reference_genome_fasta_file)
        path(reference_genes_gtf_file)
        val(params_sqanti3_qc)
        val(params_sqanti3_filter)
        val(method)
        val(output_dir)

    output:
        tuple val(sample_id), path("${method}_sqanti3/"), emit: f

    script:
        """
        mkdir -p ${method}_sqanti3/{qc,filter}/

        sqanti3_qc.py \
            --isoforms $fasta_file \
            --refGTF $reference_genes_gtf_file \
            --refFasta $reference_genome_fasta_file \
            --fasta \
            -t ${task.cpus} \
            -o ${sample_id}_${method} \
            -d ${method}_sqanti3/qc/ \
            $params_sqanti3_qc

        sqanti3_filter.py \
           $params_sqanti3_filter \
           --sqanti_class ${method}_sqanti3/qc/${sample_id}_${method}_classification.txt \
           --output ${method} \
           --dir ${method}_sqanti3/filter/
        """
}

process runSqanti3GtfMode {

    label 'sqanti3'
    tag "${sample_id}_${method}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(gtf_file)
        path(reference_genome_fasta_file)
        path(reference_genes_gtf_file)
        val(params_sqanti3_qc)
        val(params_sqanti3_filter)
        val(method)
        val(output_dir)

    output:
        tuple val(sample_id), path("${method}_sqanti3/"), emit: f

    script:
        """
        mkdir -p ${method}_sqanti3/{qc,filter}/

        sqanti3_qc.py \
            --isoforms $gtf_file \
            --refGTF $reference_genes_gtf_file \
            --refFasta $reference_genome_fasta_file \
            -t ${task.cpus} \
            -o ${sample_id}_${method} \
            -d ${method}_sqanti3/qc/ \
            $params_sqanti3_qc

        sqanti3_filter.py \
           $params_sqanti3_filter \
           --sqanti_class ${method}_sqanti3/qc/${sample_id}_${method}_classification.txt \
           --output ${method} \
           --dir ${method}_sqanti3/filter/
        """
}

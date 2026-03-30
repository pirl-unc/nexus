#!/usr/bin/env nextflow

process runEspresso {

    label 'espresso'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genes_gtf_file)
        val(params_espresso_s)
        val(params_espresso_c)
        val(params_espresso_q)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_espresso_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_espresso_outputs/
        printf "%s\t%s\n" "$bam_file" "$sample_id" >> sample.tsv

        espresso_s_script=\$(which ESPRESSO_S.pl)
        perl \$espresso_s_script \
            --list_samples sample.tsv \
            --fa $reference_genome_fasta_file \
            --anno $reference_genes_gtf_file \
            --out ${sample_id}_espresso_outputs/ \
            --num_thread ${task.cpus} \
            $params_espresso_s

        espresso_c_script=\$(which ESPRESSO_C.pl)
        perl \$espresso_c_script \
            --in ${sample_id}_espresso_outputs \
            --fa $reference_genome_fasta_file \
            --target_ID 0 \
            --num_thread ${task.cpus} \
            $params_espresso_c

        espresso_q_script=\$(which ESPRESSO_Q.pl)
        perl \$espresso_q_script \
            --list_samples ${sample_id}_espresso_outputs/sample.tsv.updated \
            --anno $reference_genes_gtf_file \
            --num_thread ${task.cpus} \
            $params_espresso_q
        """
}

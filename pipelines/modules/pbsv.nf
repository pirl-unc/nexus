#!/usr/bin/env nextflow

process runPbsv {

    label 'pbsv'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(reference_genome_fasta_file)
        val(pbsv)
        val(pbsv_discover_params)
        val(pbsv_call_params)
        val(output_dir)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_pbsv.vcf"), emit: f

    script:
        def pbsv_discover_params_ = pbsv_discover_params == true ? '' : pbsv_discover_params
        def pbsv_call_params_ = pbsv_call_params == true ? '' : pbsv_call_params

        """
        $pbsv discover $pbsv_discover_params_ $bam_file ${bam_file.baseName}_pbsv.svsig.gz
        $pbsv call \
            --num-threads ${task.cpus} \
            $pbsv_call_params_ \
            $reference_genome_fasta_file ${bam_file.baseName}_pbsv.svsig.gz ${bam_file.baseName}_pbsv.vcf
        """
}
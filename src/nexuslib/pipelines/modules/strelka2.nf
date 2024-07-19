#!/usr/bin/env nextflow

// process runStrelka2GermlineMode {
//
//     label 'strelka2'
//     tag "${sample_id}"
//     debug true
//
//    publishDir(
//         path: "${output_dir}/",
//         mode: 'copy'
//     )
//
//    input:
//         tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file)
//         val(reference_genome_fasta_file)
//         val(python2)
//         val(strelka2) // configureStrelkaGermlineWorkflow.py
//         val(strelka2_params)
//         val(output_dir)
//
//    output:
//          tuple val(sample_id), path("${sample_id}_strelka2.vcf"), emit: f
//
//    script:
//         def strelka2_params_ = strelka2_params == true ? '' : strelka2_params
//
//         """
//         mkdir -p ${sample_id}_strelka2_germline_mode/
//         ${python2} ${strelka2} \
//             --bam ${tumor_bam_file} \
//             --referenceFasta ${reference_genome_fasta_file} \
//             --runDir ${sample_id}_strelka2_germline_mode/ \
//             $strelka2_params_
//         ${python2} ${sample_id}_strelka2_germline_mode/runWorkflow.py -m local -j ${task.cpus}
//         gunzip -c ${sample_id}_strelka2_germline_mode/results/variants/variants.vcf.gz > ${sample_id}_strelka2.vcf
//         """
// }

process runStrelka2SomaticMode {

    label 'strelka2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_strelka2)
        val(output_dir)

    output:
         tuple val(sample_id), path("${sample_id}_strelka2_snvs.vcf"), path("${sample_id}_strelka2_indels.vcf"), emit: f

    script:
        """
        mkdir -p ${sample_id}_strelka2_somatic_mode/
        configureStrelkaSomaticWorkflow.py \
            --normalBam $normal_bam_file \
            --tumorBam $tumor_bam_file \
            --referenceFasta $reference_genome_fasta_file \
            --runDir ${sample_id}_strelka2_somatic_mode/ \
            $params_strelka2
        python ${sample_id}_strelka2_somatic_mode/runWorkflow.py -m local -j ${task.cpus}
        gunzip -c ${sample_id}_strelka2_somatic_mode/results/variants/somatic.snvs.vcf.gz > ${sample_id}_strelka2_snvs.vcf
        gunzip -c ${sample_id}_strelka2_somatic_mode/results/variants/somatic.indels.vcf.gz > ${sample_id}_strelka2_indels.vcf
        """
}

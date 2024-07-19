#!/usr/bin/env nextflow

process runSvaba {

    label 'svaba'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_amb_file)
        path(reference_genome_fasta_ann_file)
        path(reference_genome_fasta_bwt_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_pac_file)
        path(reference_genome_fasta_sa_file)
        val(params_svaba)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_svaba_outputs/"), emit: f

    script:
        """
        svaba run \
            --reference-genome $reference_genome_fasta_file \
            --id-string $sample_id \
            --case-bam $tumor_bam_file \
            --control-bam $normal_bam_file \
            --threads ${task.cpus} \
            $params_svaba
        mkdir -p ${sample_id}_svaba_outputs/
        mv ${sample_id}.alignments.txt.gz ${sample_id}_svaba_outputs/
        mv ${sample_id}.bps.txt.gz ${sample_id}_svaba_outputs/
        mv ${sample_id}.contigs.bam ${sample_id}_svaba_outputs/
        mv ${sample_id}.discordant.txt.gz ${sample_id}_svaba_outputs/
        mv ${sample_id}.log ${sample_id}_svaba_outputs/
        mv ${sample_id}.svaba.germline.indel.vcf ${sample_id}_svaba_outputs/
        mv ${sample_id}.svaba.germline.sv.vcf ${sample_id}_svaba_outputs/
        mv ${sample_id}.svaba.somatic.indel.vcf ${sample_id}_svaba_outputs/
        mv ${sample_id}.svaba.somatic.sv.vcf ${sample_id}_svaba_outputs/
        mv ${sample_id}.svaba.unfiltered.germline.indel.vcf ${sample_id}_svaba_outputs/
        mv ${sample_id}.svaba.unfiltered.germline.sv.vcf ${sample_id}_svaba_outputs/
        mv ${sample_id}.svaba.unfiltered.somatic.indel.vcf ${sample_id}_svaba_outputs/
        mv ${sample_id}.svaba.unfiltered.somatic.sv.vcf ${sample_id}_svaba_outputs/
        """
}

#!/usr/bin/env nextflow

process runWhatshapHaplotag {

    label 'whatshap'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(phased_vcf_file), path(phased_vcf_tbi_file)
        val(reference_genome_fasta_file)
        val(reference_genome_fasta_fai_file)
        val(params_whatshap)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_haplotagged.bam"), path("${bam_file.baseName}_haplotagged.bam.bai"), emit: f

    script:
        """
        whatshap haplotag \
            --reference $reference_genome_fasta_file \
            -o ${bam_file.baseName}_haplotagged.bam \
            $phased_vcf_file \
            $bam_file \
            $params_whatshap
        samtools index -@ ${task.cpus} -b ${bam_file.baseName}_haplotagged.bam ${bam_file.baseName}_haplotagged.bam.bai
        """
}

process runWhatshapPhase {

    label 'whatshap'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(vcf_file)
        val(reference_genome_fasta_file)
        val(reference_genome_fasta_fai_file)
        val(params_whatshap)

    output:
        tuple val(sample_id), path("${sample_id}_whatshap_phased.vcf.gz"), path("${sample_id}_whatshap_phased.vcf.gz.tbi"), emit: f

    script:
        """
        whatshap phase \
            --reference $reference_genome_fasta_file \
            --output ${sample_id}_whatshap_phased.vcf \
            $params_whatshap \
            $vcf_file $bam_file
        bgzip ${sample_id}_whatshap_phased.vcf
        tabix -p vcf ${sample_id}_whatshap_phased.vcf.gz
        """
}
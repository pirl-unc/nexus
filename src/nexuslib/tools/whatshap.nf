#!/usr/bin/env nextflow

process runWhatshapHaplotag {

    label 'whatshap'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy',
        pattern: "${bam_file.baseName}_haplotagged.bam"
    )

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy',
        pattern: "${bam_file.baseName}_haplotagged.bam.bai"
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(phased_vcf_file), path(phased_vcf_tbi_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_whatshap)
        val(output_dir)

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

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy',
        pattern: "${sample_id}_whatshap_phased.vcf.gz"
    )

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy',
        pattern: "${sample_id}_whatshap_phased.vcf.gz.tbi"
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(vcf_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_whatshap)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_whatshap_phased.vcf.gz"), path("${sample_id}_whatshap_phased.vcf.gz.tbi"), emit: f

    script:
        """
        # Ensure input VCF is bgzipped and indexed
        if [[ "${vcf_file}" == *.gz ]]; then
            vcf_gz="${vcf_file}"
        else
            bgzip -c ${vcf_file} > ${vcf_file}.gz
            vcf_gz="${vcf_file}.gz"
        fi
        tabix -p vcf \${vcf_gz}

        whatshap phase \
            --reference $reference_genome_fasta_file \
            --output ${sample_id}_whatshap_phased.vcf \
            $params_whatshap \
            \${vcf_gz} $bam_file
        bgzip ${sample_id}_whatshap_phased.vcf
        tabix -p vcf ${sample_id}_whatshap_phased.vcf.gz
        """
}

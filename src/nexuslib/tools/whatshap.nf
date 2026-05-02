#!/usr/bin/env nextflow

process runWhatshapHaplotag {

    label 'whatshap'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${bam_file.baseName}_haplotagged.bam",
        enabled: ((params.haplotag_output ?: 'bam').toString().toLowerCase() in ['bam', 'both'])
    )

    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${bam_file.baseName}_haplotagged.bam.bai",
        enabled: ((params.haplotag_output ?: 'bam').toString().toLowerCase() in ['bam', 'both'])
    )

    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${bam_file.baseName}_haplotagged_haplotag.tsv.gz",
        enabled: ((params.haplotag_output ?: 'bam').toString().toLowerCase() in ['tsv', 'both'])
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(phased_vcf_file), path(phased_vcf_tbi_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_whatshap)
        val(output_dir)
        val(subdir)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_haplotagged.bam"), path("${bam_file.baseName}_haplotagged.bam.bai"), emit: f
        path("${bam_file.baseName}_haplotagged_haplotag.tsv.gz"), emit: tsv

    script:
        // WhatsHap's --output-haplotag-list emits a 4-column TSV
        // (readname, haplotype, phaseset, chromosome). The .gz suffix triggers
        // gzipped output (per WhatsHap docs). Always generated; whether it's
        // copied to output_dir is governed by the publishDir enabled flag above.
        """
        whatshap haplotag \
            --reference $reference_genome_fasta_file \
            -o ${bam_file.baseName}_haplotagged.bam \
            --output-haplotag-list ${bam_file.baseName}_haplotagged_haplotag.tsv.gz \
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
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${sample_id}_whatshap_phased.vcf.gz"
    )

    publishDir(
        path: "${output_dir}/${sample_id}/${subdir}/",
        mode: 'copy',
        pattern: "${sample_id}_whatshap_phased.vcf.gz.tbi"
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(vcf_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        val(params_whatshap)
        val(output_dir)
        val(subdir)

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

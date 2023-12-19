#!/usr/bin/env nextflow

process runDelly2SingleSample {

    label 'delly2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'symlink'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(delly2)
        val(bcftools)
        val(reference_genome_fasta_file)
        val(exclude_regions_file)
        val(minimum_mapping_quality)
        val(output_dir)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_delly2.vcf"), path("${bam_file.baseName}_delly2.bcf"), emit: f

    script:
        """
        export OMP_NUM_THREADS=${task.cpus}
        $delly2 call \
            --exclude $exclude_regions_file \
            --outfile ${bam_file.baseName}_delly2.bcf \
            --genome $reference_genome_fasta_file \
            --map-qual $minimum_mapping_quality \
            $bam_file
        $bcftools view ${bam_file.baseName}_delly2.bcf > ${bam_file.baseName}_delly2.vcf
        """
}

process runDelly2TumorNormal {

    label 'delly2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(tumor_sample_id), val(normal_sample_id)
        val(delly2)
        val(delly2_call_params)
        val(bcftools)
        val(reference_genome_fasta_file)
        val(output_dir)

    output:
        tuple val(sample_id), path("${tumor_sample_id}_${normal_sample_id}_delly2_filtered.vcf"), path("${tumor_sample_id}_${normal_sample_id}_delly2_filtered.bcf"), emit: f

    script:
        def delly2_call_params_ = delly2_call_params == true ? '' : delly2_call_params

        """
        touch ${tumor_sample_id}_${normal_sample_id}_delly2_samples.tsv
        printf "%s\t%s\n" "$tumor_sample_id" "tumor" >> ${tumor_sample_id}_${normal_sample_id}_delly2_samples.tsv
        printf "%s\t%s\n" "$normal_sample_id" "control" >> ${tumor_sample_id}_${normal_sample_id}_delly2_samples.tsv

        export OMP_NUM_THREADS=${task.cpus}

        $delly2 call \
            --outfile ${tumor_sample_id}_${normal_sample_id}_delly2.bcf \
            --genome $reference_genome_fasta_file \
            $delly2_call_params_ \
            $tumor_bam_file $normal_bam_file

        $delly2 filter \
            -f somatic \
            -o ${tumor_sample_id}_${normal_sample_id}_delly2_filtered.bcf \
            -s ${tumor_sample_id}_${normal_sample_id}_delly2_samples.tsv \
            ${tumor_sample_id}_${normal_sample_id}_delly2.bcf

        $bcftools view ${tumor_sample_id}_${normal_sample_id}_delly2_filtered.bcf > ${tumor_sample_id}_${normal_sample_id}_delly2_filtered.vcf
        """
}
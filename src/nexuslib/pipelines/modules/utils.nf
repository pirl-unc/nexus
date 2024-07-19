#!/usr/bin/env nextflow

process copyBamFile {

    label 'copy_bam_file'
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(output_dir)

    output:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), emit: f

    script:
        """
        echo "Copying $bam_file and $bam_bai_file into $output_dir"
        mkdir -p $output_dir
        cp --remove-destination $bam_file $output_dir
        cp --remove-destination $bam_bai_file $output_dir
        """
}

process copyVcfFile {

    label 'copy_vcf_file'
    debug true

    input:
        tuple val(sample_id), path(vcf_file)
        val(output_dir)

    output:
        tuple val(sample_id), path(vcf_file), emit: f

    script:
        """
        echo "Copying $vcf_file into $output_dir"
        mkdir -p $output_dir
        cp --remove-destination $vcf_file $output_dir
        """
}

process copyIndexedVcfFile {

    label 'copy_vcf_file'
    debug true

    input:
        tuple val(sample_id), path(vcf_file), path(vcf_tbi_file)
        val(output_dir)

    output:
        tuple val(sample_id), path(vcf_file), emit: f

    script:
        """
        echo "Copying $vcf_file into $output_dir"
        mkdir -p $output_dir
        cp --remove-destination $vcf_file $output_dir
        cp --remove-destination $vcf_tbi_file $output_dir
        """
}

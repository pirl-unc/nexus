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
        echo "Copying $bam_file and $bam_bai_file into $output_dir/$sample_id/"
        mkdir -p $output_dir/$sample_id/
        cp --remove-destination $bam_file $output_dir/$sample_id/
        cp --remove-destination $bam_bai_file $output_dir/$sample_id/
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
        echo "Copying $vcf_file into $output_dir/$sample_id/"
        mkdir -p $output_dir/$sample_id/
        cp --remove-destination $vcf_file $output_dir/$sample_id/
        """
}

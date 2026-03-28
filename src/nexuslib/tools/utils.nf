#!/usr/bin/env nextflow

process decompressFile {

    label 'nexus_utils'
    debug true

    input:
        path(input_file)

    output:
        path("${decompressed_name}"), emit: f

    script:
        decompressed_name = input_file.name.endsWith('.gz') ? input_file.name[0..-4] : input_file.name
        if (input_file.name.endsWith('.gz'))
            """
            gunzip -c $input_file > ${decompressed_name}
            """
        else
            """
            echo "File is not compressed, passing through."
            """
}

process bgzipGtfFile {

    label 'nexus_utils'
    debug true

    input:
        path(gtf_file)

    output:
        tuple path("${gtf_file.baseName}.sorted.gtf.gz"), path("${gtf_file.baseName}.sorted.gtf.gz.tbi"), emit: f

    script:
        """
        sort -k1,1 -k4,4n $gtf_file | bgzip > ${gtf_file.baseName}.sorted.gtf.gz
        tabix -p gff ${gtf_file.baseName}.sorted.gtf.gz
        """
}

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
        cp $bam_file $output_dir
        sync
        cp $bam_bai_file $output_dir
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
        cp $vcf_file $output_dir
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
        cp $vcf_file $output_dir
        sync
        cp $vcf_tbi_file $output_dir
        """
}

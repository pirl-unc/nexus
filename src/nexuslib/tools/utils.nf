#!/usr/bin/env nextflow

process decompressFile {

    label 'utils_short'
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

process bgzipAndIndexVcfFile {

    label 'utils_short'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(vcf_file)

    output:
        tuple val(sample_id), path("${out_name}"), path("${out_name}.tbi"), emit: f

    script:
        out_name = vcf_file.name.endsWith('.gz') ? vcf_file.name : "${vcf_file.name}.gz"
        if (vcf_file.name.endsWith('.gz'))
            """
            tabix -f -p vcf ${vcf_file}
            """
        else
            """
            bgzip -c ${vcf_file} > ${out_name}
            tabix -f -p vcf ${out_name}
            """
}

process extractGtfFromDir {

    label 'utils_short'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(input_dir)
        val(pattern)

    output:
        tuple val(sample_id), path("extracted.gtf"), emit: f, optional: true

    script:
        // A GTF without exon features yields an empty corrected FASTA in
        // SQANTI3 GTF mode, which crashes TD2/PSAURON — leave 'extracted.gtf'
        // uncreated so the optional output drops the sample from the channel.
        """
        gtf_file=\$(find ${input_dir}/ -name "${pattern}" -type f | head -1)
        if [ -z "\$gtf_file" ]; then
            echo "WARNING: No file matching '${pattern}' found in ${input_dir}/ — skipping ${sample_id}."
        elif ! awk -F'\\t' '\$1 !~ /^#/ && \$3 == "exon" { found = 1; exit } END { exit !found }' "\$gtf_file"; then
            echo "WARNING: '\$gtf_file' has no exon features — skipping ${sample_id}."
        else
            cp "\$gtf_file" extracted.gtf
        fi
        """
}

process bgzipGtfFile {

    label 'utils_short'
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

    label 'utils_long'
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

    label 'utils_short'
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

    label 'utils_short'
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

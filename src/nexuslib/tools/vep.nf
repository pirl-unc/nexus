#!/usr/bin/env nextflow

process prepareGencodeGtfFileForVEP {

    label 'vep'
    debug true

    input:
        path(gtf_file)

    output:
        path("${gtf_file.baseName}.sorted.gtf.gz"), emit: gtf_file
        path("${gtf_file.baseName}.sorted.gtf.gz.tbi"), emit: tbi_file

    script:
        """
        grep -v "#" "${gtf_file}" \
          | sort -k1,1 -k4,4n -k5,5n -t \$'\\t' \
          | bgzip -c > ${gtf_file.baseName}.sorted.gtf.gz
        tabix -p gff ${gtf_file.baseName}.sorted.gtf.gz
        """
}

process runVEP {

    label 'vep'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(vcf_file)
        path(vep_dir)
        val(params_vep)
        val(output_dir)

    output:
        tuple val(sample_id), path("${vcf_file.baseName}_vep.tsv"), emit: f

    script:
        """
        vep \
            -i $vcf_file \
            --tab \
            --dir_cache $vep_dir \
            -o ${vcf_file.baseName}_vep.tsv \
            $params_vep
        """
}

process runVEPCustom {

    label 'vep'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(vcf_file)
        path(vep_dir)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_gtf_gz_file)
        path(reference_gtf_gz_tbi_file)
        val(reference_gtf_source)
        val(params_vep)
        val(output_dir)

    output:
        tuple val(sample_id), path("${vcf_file.baseName}_vep-custom.tsv"), emit: f

    script:
        """
        vep \
            -i $vcf_file \
            --custom file=$reference_gtf_gz_file,short_name=${reference_gtf_source},format=gtf \
            --fasta $reference_genome_fasta_file \
            --dir_cache $vep_dir \
            --tab \
            -o ${vcf_file.baseName}_vep-custom.tsv \
            $params_vep
        """
}

process runFilterVEP {

    label 'vep'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tsv_file)
        val(params_filter_vep)
        val(output_dir)

    output:
        tuple val(sample_id), path("${tsv_file.baseName}_filtered.tsv"), emit: f

    script:
        """
        filter_vep \
            -i $tsv_file \
            -o ${tsv_file.baseName}_filtered.tsv \
            $params_filter_vep
        """
}

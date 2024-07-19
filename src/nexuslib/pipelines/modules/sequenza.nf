#!/usr/bin/env nextflow

process runSequenzaUtilsBam2Seqz {

    label 'sequenza'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_wig_file)
        val(chromosomes)
        val(params_sequenza_bam2seqz)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_sequenza*seqz.gz"), path("${sample_id}_sequenza*seqz.gz.tbi"), emit: f

    script:
        """
        sequenza-utils bam2seqz \
            -t $tumor_bam_file \
            -n $normal_bam_file \
            --fasta $reference_genome_fasta_file \
            -gc $reference_genome_wig_file \
            --chromosome $chromosomes \
            --parallel ${task.cpus} \
            $params_sequenza_bam2seqz \
            -o ${sample_id}_sequenza.seqz.gz
        """
}

process mergeSequenzaSeqzFiles {

    label 'sequenza'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(seqz_gz_files), path(seqz_gz_tbi_files)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_sequenza.merged.seqz.gz"), path("${sample_id}_sequenza.merged.seqz.gz.tbi"), emit: f

    script:
        """
        zcat ${seqz_gz_files} | \
            awk '{if (NR!=1 && \$1 != "chromosome") {print \$0}}' | bgzip > \
            ${sample_id}_sequenza.merged.seqz.gz
        tabix -f -s 1 -b 2 -e 2 -S 1 ${sample_id}_sequenza.merged.seqz.gz
        """
}

process runSequenzaUtilsSeqzBinning {

    label 'sequenza'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(seqz_gz_file), path(seqz_gz_tbi_file)
        val(params_sequenza_seqzbinning)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_sequenza.small.seqz.gz"), path("${sample_id}_sequenza.small.seqz.gz.tbi"), emit: f

    script:
        """
        sequenza-utils seqz_binning \
            --seqz $seqz_gz_file \
            $params_sequenza_seqzbinning \
            -o ${sample_id}_sequenza.small.seqz.gz
        """
}

process runSequenza {

    label 'sequenza'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(small_seqz_gz_file), path(small_seqz_gz_tbi_file)
        val(chromosomes)
        val(assembly)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_sequenza/*"), emit: f

    script:
        """
        mkdir -p ${sample_id}_sequenza/
        Rscript /opt/sequenza/run_sequenza.R \
            --sample-id $sample_id \
            --small-seqz-file $small_seqz_gz_file \
            --assembly $assembly \
            --chromosomes $chromosomes \
            --output-path ${sample_id}_sequenza/
        """
}


#!/usr/bin/env nextflow

process runHLAminerShortReadDNA {

    label 'hlaminer'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/hlaminer/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(params_bwa_aln)
        val(params_bwa_sampe)
        val(params_hlaminer)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${sample_id}_HLAminer_HPRA.csv"),
              path("${sample_id}_HLAminer_HPRA.log"),
              emit: f

    script:
        """
        set -e
        bwa aln -t ${task.cpus} $params_bwa_aln \
            /opt/hlaminer/database/HLA-I_II_GEN.fasta $fastq_file_1 > aln_1.sai
        bwa aln -t ${task.cpus} $params_bwa_aln \
            /opt/hlaminer/database/HLA-I_II_GEN.fasta $fastq_file_2 > aln_2.sai
        bwa sampe $params_bwa_sampe \
            /opt/hlaminer/database/HLA-I_II_GEN.fasta \
            aln_1.sai aln_2.sai $fastq_file_1 $fastq_file_2 > aln.sam

        perl /opt/hlaminer/bin/HLAminer.pl \
            -a aln.sam \
            -h /opt/hlaminer/database/HLA-I_II_GEN.fasta \
            -p /opt/hlaminer/database/hla_nom_p.txt \
            $params_hlaminer

        mv HLAminer_HPRA.csv ${sample_id}_HLAminer_HPRA.csv
        mv HLAminer_HPRA.log ${sample_id}_HLAminer_HPRA.log
        """
}

process runHLAminerShortReadRNA {

    label 'hlaminer'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/hlaminer/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(params_bwa_aln)
        val(params_bwa_sampe)
        val(params_hlaminer)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${sample_id}_HLAminer_HPRA.csv"),
              path("${sample_id}_HLAminer_HPRA.log"),
              emit: f

    script:
        """
        set -e
        bwa aln -t ${task.cpus} $params_bwa_aln \
            /opt/hlaminer/database/HLA-I_II_CDS.fasta $fastq_file_1 > aln_1.sai
        bwa aln -t ${task.cpus} $params_bwa_aln \
            /opt/hlaminer/database/HLA-I_II_CDS.fasta $fastq_file_2 > aln_2.sai
        bwa sampe $params_bwa_sampe \
            /opt/hlaminer/database/HLA-I_II_CDS.fasta \
            aln_1.sai aln_2.sai $fastq_file_1 $fastq_file_2 > aln.sam

        perl /opt/hlaminer/bin/HLAminer.pl \
            -a aln.sam \
            -h /opt/hlaminer/database/HLA-I_II_CDS.fasta \
            -p /opt/hlaminer/database/hla_nom_p.txt \
            $params_hlaminer

        mv HLAminer_HPRA.csv ${sample_id}_HLAminer_HPRA.csv
        mv HLAminer_HPRA.log ${sample_id}_HLAminer_HPRA.log
        """
}

process runHLAminerLongReadDNA {

    label 'hlaminer'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/hlaminer/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_minimap2)
        val(params_hlaminer)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${sample_id}_HLAminer_HPRA.csv"),
              path("${sample_id}_HLAminer_HPRA.log"),
              emit: f

    script:
        """
        set -o pipefail
        minimap2 \
            -t ${task.cpus} \
            --MD \
            $params_minimap2 \
            /opt/hlaminer/database/HLA-I_II_GEN.fasta \
            $fastq_file \
            | /opt/hlaminer/bin/HLAminer.pl \
                -a stream \
                -h /opt/hlaminer/database/HLA-I_II_GEN.fasta \
                -p /opt/hlaminer/database/hla_nom_p.txt \
                $params_hlaminer

        mv HLAminer_HPRA.csv ${sample_id}_HLAminer_HPRA.csv
        mv HLAminer_HPRA.log ${sample_id}_HLAminer_HPRA.log
        """
}

process runHLAminerLongReadRNA {

    label 'hlaminer'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/hlaminer/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_minimap2)
        val(params_hlaminer)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${sample_id}_HLAminer_HPRA.csv"),
              path("${sample_id}_HLAminer_HPRA.log"),
              emit: f

    script:
        """
        set -o pipefail
        minimap2 \
            -t ${task.cpus} \
            --MD \
            $params_minimap2 \
            /opt/hlaminer/database/HLA-I_II_CDS.fasta \
            $fastq_file \
            | /opt/hlaminer/bin/HLAminer.pl \
                -a stream \
                -h /opt/hlaminer/database/HLA-I_II_CDS.fasta \
                -p /opt/hlaminer/database/hla_nom_p.txt \
                $params_hlaminer

        mv HLAminer_HPRA.csv ${sample_id}_HLAminer_HPRA.csv
        mv HLAminer_HPRA.log ${sample_id}_HLAminer_HPRA.log
        """
}
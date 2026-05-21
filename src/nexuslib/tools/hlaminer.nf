#!/usr/bin/env nextflow

process runHLAminer {

    label 'hlaminer'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/hlaminer/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(params_hlaminer)
        val(output_dir)

    output:
        tuple val(sample_id),
              path("${sample_id}_HLAminer_HPRA.csv"),
              path("${sample_id}_HLAminer_HPRA.log"),
              emit: f

    script:
        """
        samtools view $bam_file \
            | perl /opt/hlaminer/bin/HLAminer.pl \
                -a /dev/stdin \
                -h /opt/hlaminer/database/HLA-I_II_CDS.fasta.gz \
                -p /opt/hlaminer/database/hla_nom_p.txt \
                $params_hlaminer

        # Prefix HLAminer's hardcoded output filenames with sample_id so the
        # emitted tuple is keyed per-sample and publishDir doesn't collide
        # if multiple samples ever land in the same target dir.
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
        # Align long-read RNA reads against HLAminer's HLA CDS reference
        # and pipe the SAM straight into HLAminer.pl. No intermediate BAM is
        # written; mirrors HLAminer's official wrapper-script convention
        # (HPTASRrnaseq_classI-II.sh etc.).
        #
        # The HLA CDS reference contains transcript-style sequences (no
        # introns), so a non-splice long-read preset (e.g. -ax map-hifi for
        # PacBio HiFi, -ax map-ont for Nanopore) is the right choice in
        # params_minimap2.
        #
        # set -o pipefail so a minimap2 failure isn't silently masked by
        # HLAminer.pl's zero exit (HLAminer.pl exits 0 even on empty input).
        set -o pipefail
        minimap2 \
            -t ${task.cpus} \
            $params_minimap2 \
            /opt/hlaminer/database/HLA-I_II_CDS.fasta.gz \
            $fastq_file \
            | /opt/hlaminer/bin/HLAminer.pl \
                -a stream \
                -h /opt/hlaminer/database/HLA-I_II_CDS.fasta.gz \
                -p /opt/hlaminer/database/hla_nom_p.txt \
                $params_hlaminer

        # Prefix HLAminer's hardcoded output filenames with sample_id so the
        # emitted tuple is keyed per-sample and publishDir doesn't collide
        # if multiple samples ever land in the same target dir.
        mv HLAminer_HPRA.csv ${sample_id}_HLAminer_HPRA.csv
        mv HLAminer_HPRA.log ${sample_id}_HLAminer_HPRA.log
        """
}


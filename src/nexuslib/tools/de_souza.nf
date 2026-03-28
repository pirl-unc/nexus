#!/usr/bin/env nextflow

process runFlagCorrection {

    label 'de_souza'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(splitncigarreads_bam_file), path(splitncigarreads_bam_bai_file)

    output:
        tuple val(sample_id), path("${splitncigarreads_bam_file.baseName}_flagcorrected.bam"), emit: f

    script:
        """
        Rscript /opt/lrRNAseqVariantCalling/flagCorrection.r \
            $bam_file \
            $splitncigarreads_bam_file \
            ${splitncigarreads_bam_file.baseName}_flagcorrected.bam \
            ${task.cpus}
        """
}

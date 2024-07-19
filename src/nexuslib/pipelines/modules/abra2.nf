#!/usr/bin/env nextflow

process runAbra2 {

    label 'abra2'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(abra2_targets_bed_file)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_realigned.bam"), path("${bam_file.baseName}_realigned.bam.bai"), emit: f

    script:
        """
        mkdir -p abra2_temp/
        java -Xmx${task.java_max_mem.toGiga()}G -jar /opt/abra2/abra2.jar \
            --in $bam_file \
            --out ${bam_file.baseName}_realigned.bam \
            --ref $reference_genome_fasta_file \
            --targets $abra2_targets_bed_file \
            --threads ${task.cpus} \
            --tmpdir abra2_temp/
        samtools index -b ${bam_file.baseName}_realigned.bam ${bam_file.baseName}_realigned.bam.bai
        """
}
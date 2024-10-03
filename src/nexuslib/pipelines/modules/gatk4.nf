#!/usr/bin/env nextflow

process runGatk4BaseRecalibrator {

    label 'gatk4_baserecalibrator'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), val(chromosome)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_gzi_file)
        path(reference_genome_fasta_dict_file)
        val(known_sites_files_args)

    output:
         tuple val(sample_id), path(bam_file), val(chromosome), path("${sample_id}_recalibration_${chromosome}.data.table"), emit: f

    script:
        """
        gatk --java-options -Xmx${task.java_max_mem.toGiga()}G BaseRecalibrator \
            -I $bam_file \
            -L $chromosome \
            -R $reference_genome_fasta_file \
            -O ${sample_id}_recalibration_${chromosome}.data.table \
            $known_sites_files_args
        """
}

process runGatk4GatherBQSRReports {

    label 'gatk4_gatherbqsrreports'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(data_table_files)

    output:
        tuple val(sample_id), path("${sample_id}_recalibration_merged.data.table"), emit: f

    script:
        """
        data_table_files_str="${data_table_files}"
        IFS=' ' array=(\${data_table_files_str})
        input_data_table_files=""
        for element in "\${array[@]}"; do
            input_data_table_files+="-I \${element} "
        done
        gatk --java-options -Xmx${task.java_max_mem.toGiga()}G GatherBQSRReports \
            \${input_data_table_files} \
            -O ${sample_id}_recalibration_merged.data.table
        """
}

process runGatk4ApplyBQSR {

    label 'gatk4_applybqsr'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(data_table_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_gzi_file)
        path(reference_genome_fasta_dict_file)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_recalibrated.bam"), path("${bam_file.baseName}_recalibrated.bam.bai"), emit: f

    script:
        """
        gatk --java-options -Xmx${task.java_max_mem.toGiga()}G ApplyBQSR \
            -I ${bam_file} \
            -R ${reference_genome_fasta_file} \
            --bqsr-recal-file ${data_table_file} \
            -O ${bam_file.baseName}_recalibrated.bam
        samtools index -b ${bam_file.baseName}_recalibrated.bam ${bam_file.baseName}_recalibrated.bam.bai
        """
}

process runGatk4LearnReadOrientationModel {

    label 'gatk4_learnreadorientationmodel'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), val(chromosome), path(f1r2_file)

    output:
        tuple val(sample_id), val(chromosome), path("${f1r2_file.baseName}_gatk4_read-orientation-model.tar.gz"), emit: f

    script:
        """
        gatk --java-options -Xmx${task.java_max_mem.toGiga()}G LearnReadOrientationModel \
            -I $f1r2_file \
            -O ${f1r2_file.baseName}_gatk4_read-orientation-model.tar.gz
        """
}

process runGatk4GetPileupSummaries {

    label 'gatk4_getpileupsummaries'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(params_gatk4getpileupsummaries)

    output:
        tuple val(sample_id), path("${sample_id}_gatk4_pileup-summaries.table"), emit: f

    script:
        """
        gatk --java-options -Xmx${task.java_max_mem.toGiga()}G GetPileupSummaries \
            -I $bam_file \
            -O ${sample_id}_gatk4_pileup-summaries.table \
            $params_gatk4getpileupsummaries
        """
}

process runGatk4CalculateContamination {

    label 'gatk4_calculatecontamination'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(pileup_file)

    output:
        tuple val(sample_id), path("${sample_id}_gatk4_segments.table"), path("${sample_id}_gatk4_contamination.table"), emit: f

    script:
        """
        gatk --java-options -Xmx${task.java_max_mem.toGiga()}G CalculateContamination \
            -I $pileup_file \
            -tumor-segmentation ${sample_id}_gatk4_segments.table \
            -O ${sample_id}_gatk4_contamination.table
        """
}

process runGatk4FilterMutect2Calls {

    // Use this process when both segmentation file and contamination file were generated

    label 'gatk4_filtermutectcalls'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), val(chromosome), path(vcf_file), path(vcf_idx_file), path(vcf_stats_file), path(segmentation_file), path(contamination_file), path(read_orientation_model_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_dict_file)

    output:
        tuple val(sample_id), path("${sample_id}_gatk4-mutect2_${chromosome}_filtered.vcf"), emit: f

    script:
        """
        gatk --java-options -Xmx${task.java_max_mem.toGiga()}G FilterMutectCalls \
            -V $vcf_file \
            -R $reference_genome_fasta_file \
            --tumor-segmentation $segmentation_file \
            --contamination-table $contamination_file \
            --ob-priors $read_orientation_model_file \
            -O ${sample_id}_gatk4-mutect2_${chromosome}_filtered.vcf
        """
}

process runGatk4FilterMutect2CallsNonHumanSample {

    // Use this process when neither segmentation file nor contamination file was generated (i.e. BAM file is for a non-human sample)

    label 'gatk4_filtermutectcalls'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), val(chromosome), path(vcf_file), path(vcf_idx_file), path(vcf_stats_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_dict_file)

    output:
        tuple val(sample_id), path("${sample_id}_gatk4-mutect2_${chromosome}_filtered.vcf"), emit: f

    script:
        """
        gatk --java-options -Xmx${task.java_max_mem.toGiga()}G FilterMutectCalls \
            -V $vcf_file \
            -R $reference_genome_fasta_file \
            -O ${sample_id}_gatk4-mutect2_${chromosome}_filtered.vcf
        """
}

process runGatk4Mutect2TumorNormal {

    label 'gatk4_mutect2'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(normal_sample_id), val(chromosome)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_dict_file)
        val(params_gatk4mutect2)

    output:
        tuple val(sample_id), val(chromosome), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf.idx"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf.stats"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered_f1r2.tar.gz"), emit: f

    script:
        """
        gatk --java-options -Xmx${task.java_max_mem.toGiga()}G Mutect2 \
            -R $reference_genome_fasta_file \
            -I $tumor_bam_file \
            -I $normal_bam_file \
            -normal $normal_sample_id \
            --native-pair-hmm-threads ${task.hmm_threads} \
            --intervals $chromosome \
            --f1r2-tar-gz ${sample_id}_gatk4-mutect2_${chromosome}_unfiltered_f1r2.tar.gz \
            -O ${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf \
            $params_gatk4mutect2
        """
}

process runGatk4Mutect2TumorOnly {

    label 'gatk4_mutect2'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), val(chromosome)
        path(reference_genome_fasta_file)
        val(params_gatk4mutect2)

    output:
        tuple val(sample_id), val(chromosome), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf.idx"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf.stats"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered_f1r2.tar.gz"), emit: f

    script:
        """
        gatk --java-options -Xmx${task.java_max_mem.toGiga()}G Mutect2 \
            -R $reference_genome_fasta_file \
            -I $tumor_bam_file \
            --native-pair-hmm-threads ${task.hmm_threads} \
            --intervals $chromosome \
            --f1r2-tar-gz ${sample_id}_gatk4-mutect2_${chromosome}_unfiltered_f1r2.tar.gz \
            -O ${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf \
            $params_gatk4mutect2
        """
}

process runGatk4SplitNCigarReads {

    label 'gatk4_splitncigarreads'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_dict_file)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_gatk4_splitncigarreads.bam"), emit: f

    script:
        """
        gatk --java-options "-Xmx${task.java_max_mem.toGiga()}G -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}" SplitNCigarReads \
            -R ${reference_genome_fasta_file} \
            -I ${bam_file} \
            -O ${bam_file.baseName}_gatk4_splitncigarreads.bam
        """
}

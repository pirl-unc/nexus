#!/usr/bin/env nextflow

process runGatk4BaseRecalibrator {

    label 'gatk4_base_recalibrator'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), val(chromosome)
        val(reference_genome_fasta_file)
        val(gatk4)
        val(gatk4_baserecalibrator_params)

    output:
         tuple val(sample_id), path(bam_file), val(chromosome), path("${sample_id}_recalibration_${chromosome}.data.table"), emit: f

    script:
        def gatk4_baserecalibrator_params_ = gatk4_baserecalibrator_params == true ? '' : gatk4_baserecalibrator_params

        """
        $gatk4 BaseRecalibrator \
            -I $bam_file \
            -L $chromosome \
            -R $reference_genome_fasta_file \
            -O ${sample_id}_recalibration_${chromosome}.data.table \
            $gatk4_baserecalibrator_params_
        """
}

process runGatk4GatherBQSRReports {

    label 'gatk4_gather_bqsr_reports'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(data_table_files)
        val(gatk4)

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
        CMD="${gatk4} --java-options -Xmx${task.java_max_mem.toGiga()}G GatherBQSRReports \${input_data_table_files} -O ${sample_id}_recalibration_merged.data.table"
        eval \${CMD}
        """
}

process runGatk4ApplyBQSR {

    label 'gatk4_apply_bqsr'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(data_table_file)
        val(reference_genome_fasta_file)
        val(gatk4)
        val(samtools)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_recalibrated.bam"), path("${bam_file.baseName}_recalibrated.bam.bai"), emit: f

    script:
        """
        $gatk4 --java-options -Xmx${task.java_max_mem.toGiga()}G ApplyBQSR \
            -I ${bam_file} \
            -R ${reference_genome_fasta_file} \
            --bqsr-recal-file ${data_table_file} \
            -O ${bam_file.baseName}_recalibrated.bam
        $samtools index -b ${bam_file.baseName}_recalibrated.bam ${bam_file.baseName}_recalibrated.bam.bai
        """
}

process runGatk4LearnReadOrientationModel {

    label 'gatk4_learn_read_orientation_model'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), val(chromosome), path(f1r2_file)
        val(gatk4)

    output:
        tuple val(sample_id), val(chromosome), path("${f1r2_file.baseName}_gatk4_read-orientation-model.tar.gz"), emit: f

    script:
        """
        $gatk4 --java-options -Xmx${task.java_max_mem.toGiga()}G LearnReadOrientationModel \
            -I $f1r2_file \
            -O ${f1r2_file.baseName}_gatk4_read-orientation-model.tar.gz
        """
}

process runGatk4GetPileupSummaries {

    label 'gatk4_get_pileup_summaries'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(gatk4)
        val(gatk4_getpileupsummaries_params)

    output:
        tuple val(sample_id), path("${sample_id}_gatk4_pileup-summaries.table"), emit: f

    script:
        def gatk4_getpileupsummaries_params_ = gatk4_getpileupsummaries_params == true ? '' : gatk4_getpileupsummaries_params

        """
        $gatk4 --java-options -Xmx${task.java_max_mem.toGiga()}G GetPileupSummaries \
            -I $bam_file \
            -O ${sample_id}_gatk4_pileup-summaries.table \
            $gatk4_getpileupsummaries_params_
        """
}

process runGatk4CalculateContamination {

    label 'gatk4_calculate_contamination'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(pileup_file)
        val(gatk4)

    output:
        tuple val(sample_id), path("${sample_id}_gatk4_segments.table"), path("${sample_id}_gatk4_contamination.table"), emit: f

    script:
        """
        $gatk4 --java-options -Xmx${task.java_max_mem.toGiga()}G CalculateContamination \
            -I $pileup_file \
            -tumor-segmentation ${sample_id}_gatk4_segments.table \
            -O ${sample_id}_gatk4_contamination.table
        """
}

process runGatk4FilterMutect2Calls {

    // Use this process when both segmentation file and contamination file were generated

    label 'gatk4_filter_mutect_calls'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), val(chromosome), path(vcf_file), path(vcf_idx_file), path(vcf_stats_file), path(segmentation_file), path(contamination_file), path(read_orientation_model_file)
        val(reference_genome_fasta_file)
        val(gatk4)

    output:
        tuple val(sample_id), path("${sample_id}_gatk4-mutect2_${chromosome}_filtered.vcf"), emit: f

    script:
        """
        $gatk4 --java-options -Xmx${task.java_max_mem.toGiga()}G FilterMutectCalls \
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

    label 'gatk4_filter_mutect_calls'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), val(chromosome), path(vcf_file), path(vcf_idx_file), path(vcf_stats_file)
        val(gatk4)
        val(reference_genome_fasta_file)

    output:
        tuple val(sample_id), path("${sample_id}_gatk4-mutect2_${chromosome}_filtered.vcf"), emit: f

    script:
        """
        $gatk4 --java-options -Xmx${task.java_max_mem.toGiga()}G FilterMutectCalls \
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
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(tumor_sample_id), val(normal_sample_id), val(chromosome)
        val(reference_genome_fasta_file)
        val(gatk4)
        val(gat4k_mutect2_params)

    output:
        tuple val(sample_id), val(chromosome), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf.idx"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf.stats"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered_f1r2.tar.gz"), emit: f

    script:
        def gat4k_mutect2_params_ = gat4k_mutect2_params == true ? '' : gat4k_mutect2_params

        """
        $gatk4 --java-options -Xmx${task.java_max_mem.toGiga()}G Mutect2 \
            -R $reference_genome_fasta_file \
            -I $tumor_bam_file \
            -I $normal_bam_file \
            -normal $normal_sample_id \
            --native-pair-hmm-threads ${task.hmm_threads} \
            --intervals $chromosome \
            --f1r2-tar-gz ${sample_id}_gatk4-mutect2_${chromosome}_unfiltered_f1r2.tar.gz \
            -O ${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf \
            $gat4k_mutect2_params_
        """
}

process runGatk4Mutect2TumorOnly {

    label 'gatk4_mutect2'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), val(chromosome)
        val(reference_genome_fasta_file)
        val(gatk4)
        val(gat4k_mutect2_params)

    output:
        tuple val(sample_id), val(chromosome), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf.idx"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf.stats"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered_f1r2.tar.gz"), emit: f

    script:
        def gat4k_mutect2_params_ = gat4k_mutect2_params == true ? '' : gat4k_mutect2_params

        """
        $gatk4 --java-options -Xmx${task.java_max_mem.toGiga()}G Mutect2 \
            -R $reference_genome_fasta_file \
            -I $tumor_bam_file \
            --native-pair-hmm-threads ${task.hmm_threads} \
            --intervals $chromosome \
            --f1r2-tar-gz ${sample_id}_gatk4-mutect2_${chromosome}_unfiltered_f1r2.tar.gz \
            -O ${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf \
            $gat4k_mutect2_params_
        """
}

process runGatk4Mutect2TumorNormalNonHumanSample {

    label 'gatk4_mutect2'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(tumor_sample_id), val(normal_sample_id), val(chromosome)
        val(gatk4)
        val(reference_genome_fasta_file)

    output:
        tuple val(sample_id), val(chromosome), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf.idx"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf.stats"), path("${sample_id}_gatk4-mutect2_${chromosome}_unfiltered_f1r2.tar.gz"), emit: f

    script:
        """
        $gatk4 --java-options -Xmx${task.java_max_mem.toGiga()}G Mutect2 \
            -R $reference_genome_fasta_file \
            -I $tumor_bam_file \
            -I $normal_bam_file \
            -normal $normal_sample_id \
            --native-pair-hmm-threads ${task.hmm_threads} \
            --intervals $chromosome \
            --f1r2-tar-gz ${sample_id}_gatk4-mutect2_${chromosome}_unfiltered_f1r2.tar.gz \
            -O ${sample_id}_gatk4-mutect2_${chromosome}_unfiltered.vcf
        """
}
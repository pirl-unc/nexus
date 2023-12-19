#!/usr/bin/env nextflow

process runExactoConvert {

    label 'exacto_convert'
    tag "${sample_id}"
    debug true

   publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(vcf_file)
        val(exacto)
        val(exacto_convert_params)
        val(variant_calling_method)
        val(sequencing_platform)
        val(output_dir)

    output:
        tuple val(sample_id), path("${vcf_file.baseName}.tsv"), emit: f

    script:
        def exacto_convert_params_ = exacto_convert_params == true ? '' : exacto_convert_params

        """
        $exacto convert \
            --vcf-file $vcf_file \
            --variant-calling-method $variant_calling_method \
            --source-id $sample_id \
            --output-tsv-file ${vcf_file.baseName}.tsv \
            $exacto_convert_params_
        """
}

process runExactoRefineSingleSample {

    label 'exacto_refine'
    tag "${sample_id}"
    debug true

   publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(vcf_file), val(tumor_sample_id), val(normal_sample_id)
        val(exacto)
        val(exacto_refine_params)
        val(variant_calling_method)
        val(output_dir)

    output:
        tuple val(sample_id), path("${vcf_file.baseName}_refined.tsv"), emit: f

    script:
        def exacto_refine_params_ = exacto_refine_params == true ? '' : exacto_refine_params
        def tumor_sample_id_ = tumor_sample_id == true ? '' : tumor_sample_id
        def normal_sample_id_ = normal_sample_id == true ? '' : normal_sample_id

        """
        $exacto refine \
            --vcf_file $vcf_file \
            --variant_calling_method $variant_calling_method \
            --num_processes ${task.cpus} \
            --sample_id $sample_id \
            --tumor_sample_id $tumor_sample_id_ \
            --normal_sample_id $normal_sample_id_ \
            --output_tsv_file ${vcf_file.baseName}_refined.tsv \
            $exacto_refine_params_
        """
}

process runExactoRefineSomaticSnvIndel {

    label 'exacto_refine'
    tag "${sample_id}"
    debug true

   publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_vcf_file), path(normal_tsv_file), val(tumor_sample_id), val(normal_sample_id)
        val(exacto)
        val(exacto_refine_params)
        val(variant_calling_method)
        val(output_dir)

    output:
        tuple val(sample_id), path("${tumor_vcf_file.baseName}_refined.tsv"), emit: f

    script:
        def exacto_refine_params_ = exacto_refine_params == true ? '' : exacto_refine_params

        """
        $exacto refine \
            --vcf_file $tumor_vcf_file \
            --variant_calling_method $variant_calling_method \
            --num_processes ${task.cpus} \
            --sample_id $sample_id \
            --tumor_sample_id $tumor_sample_id \
            --exclude_snv_indel_tsv_files $normal_tsv_file \
            --output_tsv_file ${tumor_vcf_file.baseName}_refined.tsv \
            $exacto_refine_params_
        """
}

process runExactoRefineSomaticSv {

    label 'exacto_refine'
    tag "${sample_id}"
    debug true

   publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_vcf_file), path(normal_tsv_file), val(tumor_sample_id), val(normal_sample_id)
        val(exacto)
        val(exacto_refine_params)
        val(variant_calling_method)
        val(output_dir)

    output:
        tuple val(sample_id), path("${tumor_vcf_file.baseName}_refined.tsv"), emit: f

    script:
        def exacto_refine_params_ = exacto_refine_params == true ? '' : exacto_refine_params

        """
        $exacto refine \
            --vcf_file $tumor_vcf_file \
            --variant_calling_method $variant_calling_method \
            --num_processes ${task.cpus} \
            --sample_id $sample_id \
            --tumor_sample_id $tumor_sample_id \
            --exclude_sv_tsv_files $normal_tsv_file \
            --output_tsv_file ${tumor_vcf_file.baseName}_refined.tsv \
            $exacto_refine_params_
        """
}

process runExactoAnnotate {

    label 'exacto_annotate'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tsv_file)
        val(exacto)
        val(annotation_source)
        val(variant_class)
        val(ensembl_release)
        val(ensembl_species)
        val(gencode_gtf_file)
        val(perl_path)
        val(annovar_path)
        val(annovar_humandb_path)
        val(annovar_genome_assembly)
        val(annovar_protocol)
        val(annovar_operation)
        val(output_dir)

    output:
        tuple val(sample_id), path("${tsv_file.baseName}_${annotation_source}_annotated.tsv"), emit: f

    script:
        def ensembl_release_ = ensembl_release != "" ? "--ensembl_release $ensembl_release" : ''
        def ensembl_species_ = ensembl_species != "" ? "--ensembl_species $ensembl_species" : ''
        def gencode_gtf_file_ = gencode_gtf_file != "" ? "--gencode_gtf_file $gencode_gtf_file" : ''
        def perl_path_ = perl_path != "" ? "--perl_path $perl_path" : ''
        def annovar_path_ = annovar_path != "" ? "--annovar_path $annovar_path" : ''
        def annovar_humandb_path_ = annovar_humandb_path != "" ? "--annovar_humandb_path $annovar_humandb_path" : ''
        def annovar_genome_assembly_ = annovar_genome_assembly != "" ? "--annovar_genome_assembly $annovar_genome_assembly" : ''
        def annovar_protocol_ = annovar_protocol != "" ? "--annovar_protocol '$annovar_protocol'" : ''
        def annovar_operation_ = annovar_operation != "" ? "--annovar_operation '$annovar_operation'" : ''

        """
        $exacto annotate \
            --annotation_source $annotation_source \
            --variant_class $variant_class \
            --tsv_file $tsv_file \
            --output_tsv_file ${tsv_file.baseName}_${annotation_source}_annotated.tsv \
            --output_avinput_file ${tsv_file.baseName}.avinput \
            $ensembl_release_ \
            $ensembl_species_ \
            $gencode_gtf_file_ \
            $perl_path_ \
            $annovar_path_ \
            $annovar_humandb_path_ \
            $annovar_genome_assembly_ \
            $annovar_protocol_ \
            $annovar_operation_
        """
}

process runExactoMergeVariants {

    label 'exacto_merge_variants'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tsv_files), val(callset_sample_id)
        val(exacto)
        val(variant_class)
        val(max_clustering_distance)
        val(output_dir)

    output:
        tuple val(sample_id), path("${callset_sample_id}_${variant_class}_merged_deduped.tsv"), emit: f

    script:
        """
        tsv_files_str="${tsv_files}"
        IFS=' ' array=(\${tsv_files_str})
        input_tsv_files_arr=""
        for element in "\${array[@]}"; do
            input_tsv_files_arr+=(\${element})
        done

        $exacto merge-variants \
            --variant_class $variant_class \
            --tsv_files \${input_tsv_files_arr[@]} \
            --output_merged_tsv_file ${callset_sample_id}_${variant_class}_merged.tsv \
            --output_merged_deduped_tsv_file ${callset_sample_id}_${variant_class}_merged_deduped.tsv \
            --max_clustering_distance $max_clustering_distance
        """
}

process runExactoMergeAnnotations {

    label 'exacto_merge_annotations'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tsv_files)
        val(exacto)
        val(output_file_suffix)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_${output_file_suffix}.tsv"), emit: f

    script:
        """
        tsv_files_str="${tsv_files}"
        IFS=' ' array=(\${tsv_files_str})
        input_tsv_files_arr=""
        for element in "\${array[@]}"; do
            input_tsv_files_arr+=(\${element})
        done

        $exacto merge-annotations \
            --tsv_files \${input_tsv_files_arr[@]} \
            --output_merged_tsv_file ${sample_id}_${output_file_suffix}.tsv
        """
}

process runExactoCallVariants {

    label 'exacto_call_variants'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(exacto)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_${output_file_suffix}.tsv"), emit: f

    script:
        """
        exacto call-variants \


        tsv_files_str="${tsv_files}"
        IFS=' ' array=(\${tsv_files_str})
        input_tsv_files_arr=""
        for element in "\${array[@]}"; do
            input_tsv_files_arr+=(\${element})
        done

        $exacto merge-annotations \
            --tsv_files \${input_tsv_files_arr[@]} \
            --output_merged_tsv_file ${sample_id}_${output_file_suffix}.tsv
        """
}
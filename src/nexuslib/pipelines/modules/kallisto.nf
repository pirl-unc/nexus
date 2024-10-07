#!/usr/bin/env nextflow

process runKallistoQuantTccLongReads {

    label 'kallisto'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        path(kallisto_index_file)
        path(gtf_file)
        path(t2g_file)
        val(params_kallisto_bus)
        val(params_bustools_sort)
        val(params_bustools_count)
        val(params_kallisto_quanttcc)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_kallisto_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_kallisto_outputs/

        kallisto bus \
            -t ${task.cpus} \
            --long \
            -i $kallisto_index_file \
            -o ${sample_id}_kallisto_outputs/ \
            $params_kallisto_bus \
            $fastq_file

        bustools sort \
            -t ${task.cpus} \
            -o ${sample_id}_kallisto_outputs/sorted.bus \
            -m ${task.bustools_memory.toGiga()}G \
            $params_bustools_sort \
            ${sample_id}_kallisto_outputs/output.bus

        bustools count \
            -t ${sample_id}_kallisto_outputs/transcripts.txt \
            -e ${sample_id}_kallisto_outputs/matrix.ec \
            -o ${sample_id}_kallisto_outputs/count \
            -g $t2g_file \
            $params_bustools_count \
            ${sample_id}_kallisto_outputs/sorted.bus

        kallisto quant-tcc \
            -t ${task.cpus} \
            --long \
            -f ${sample_id}_kallisto_outputs/flens.txt \
            -i $kallisto_index_file \
            -e ${sample_id}_kallisto_outputs/count.ec.txt \
            -o ${sample_id}_kallisto_outputs/ \
            --gtf $gtf_file \
            $params_kallisto_quanttcc \
            ${sample_id}_kallisto_outputs/count.mtx
        """
}

process runKallistoQuantSingleEndReads {

    label 'kallisto'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        path(kallisto_index_file)
        val(params_kallisto_quant_fragment_length)
        val(params_kallisto_quant_sd)
        val(params_kallisto_quant)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_kallisto_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_kallisto_outputs/
        kallisto quant \
            --index=$kallisto_index_file \
            --output-dir=${sample_id}_kallisto_outputs/ \
            --fragment-length=$params_kallisto_quant_fragment_length \
            --sd=$params_kallisto_quant_sd \
            --single \
            --threads=${task.cpus} \
            $params_kallisto_quant \
            $fastq_file
        """
}

process runKallistoQuantPairedEndReads {

    label 'kallisto'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        path(kallisto_index_file)
        val(params_kallisto_quant)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_kallisto_outputs/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_kallisto_outputs/
        kallisto quant \
            --index=$kallisto_index_file \
            --output-dir=${sample_id}_kallisto_outputs/ \
            --threads=${task.cpus} \
            $params_kallisto_quant \
            $fastq_file_1 $fastq_file_2
        """
}

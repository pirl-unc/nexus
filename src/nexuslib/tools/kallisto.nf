#!/usr/bin/env nextflow

process runKallistoIndex {

    label 'kallisto_index'
    debug true

    input:
        path(reference_transcripts_fasta_file)
        val(params_kallisto_index)

    output:
        path("${reference_transcripts_fasta_file.baseName}_kallisto.index"), emit: index_file

    script:
        """
        kallisto index \
            --index ${reference_transcripts_fasta_file.baseName}_kallisto.index \
            --threads ${task.cpus} \
            $params_kallisto_index \
            $reference_transcripts_fasta_file
        """
}

process runKallistoT2G {

    label 'kallisto_index'
    debug true

    input:
        path(reference_genes_gtf_file)

    output:
        path("${reference_genes_gtf_file.baseName}_kallisto.t2g"), emit: t2g_file

    script:
        """
        t2g.py $reference_genes_gtf_file > ${reference_genes_gtf_file.baseName}_kallisto.t2g
        """
}

process runKallistoQuantTccLongReads {

    label 'kallisto'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file)
        path(kallisto_index_file)
        path(reference_genes_gtf_file)
        path(reference_genes_t2g_file)
        val(params_kallisto_bus)
        val(params_bustools_sort)
        val(params_bustools_count)
        val(params_kallisto_quanttcc)
        val(output_dir)

    output:
        tuple val(sample_id), path("kallisto/"), emit: f

    script:
        """
        mkdir -p kallisto/

        kallisto bus \
            -t ${task.cpus} \
            --long \
            -i $kallisto_index_file \
            -o kallisto/ \
            $params_kallisto_bus \
            $fastq_file

        bustools sort \
            -t ${task.cpus} \
            -o kallisto/sorted.bus \
            -m ${task.bustools_memory.toGiga()}G \
            $params_bustools_sort \
            kallisto/output.bus

        bustools count \
            -t kallisto/transcripts.txt \
            -e kallisto/matrix.ec \
            -o kallisto/count \
            -g $reference_genes_t2g_file \
            $params_bustools_count \
            kallisto/sorted.bus

        kallisto quant-tcc \
            -t ${task.cpus} \
            --long \
            -f kallisto/flens.txt \
            -i $kallisto_index_file \
            -e kallisto/count.ec.txt \
            -o kallisto/ \
            --gtf $reference_genes_gtf_file \
            $params_kallisto_quanttcc \
            kallisto/count.mtx
        """
}

process runKallistoQuantSingleEndReads {

    label 'kallisto'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
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
        tuple val(sample_id), path("kallisto/"), emit: f

    script:
        """
        mkdir -p kallisto/
        kallisto quant \
            --index=$kallisto_index_file \
            --output-dir=kallisto/ \
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
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        path(kallisto_index_file)
        val(params_kallisto_quant)
        val(output_dir)

    output:
        tuple val(sample_id), path("kallisto/"), emit: f

    script:
        """
        mkdir -p kallisto/
        kallisto quant \
            --index=$kallisto_index_file \
            --output-dir=kallisto/ \
            --threads=${task.cpus} \
            $params_kallisto_quant \
            $fastq_file_1 $fastq_file_2
        """
}

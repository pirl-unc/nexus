#!/usr/bin/env nextflow

// process runSalmonAlignmentMode {
//
//     label 'salmon'
//     tag "${sample_id}"
//     debug true
//
//     publishDir(
//         path: "${output_dir}/${sample_id}/",
//         mode: 'copy'
//     )
//
//     input:
//         tuple val(sample_id), path(bam_file), path(bam_bai_file)
//         val(salmon)
//         val(salmon_params)
//         val(output_dir)
//
//     output:
//         tuple val(sample_id), path("salmon/"), emit: f
//
//     script:
//         """
//         mkdir -p salmon/
//         $salmon quant \
//             --alignments $bam_file \
//             --output salmon/ \
//             --threads ${task.cpus} \
//             $salmon_params
//         """
// }

process runSalmonPairedEndMappingMode {

    label 'salmon'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        path(reference_transcripts_fasta_file)
        path(gtf_file)
        val(params_salmon_index)
        val(params_salmon_quant)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_salmon_output/"), emit: f

    script:
        """
        mkdir -p ${sample_id}_salmon_index/
        salmon index \
            --transcripts $reference_transcripts_fasta_file \
            --index ${sample_id}_salmon_index/ \
            --threads ${task.cpus} \
            $params_salmon_index

        mkdir -p ${sample_id}_salmon_output/
        salmon quant \
            $params_salmon_quant \
            --index ${sample_id}_salmon_index/ \
            --mates1 $fastq_file_1 \
            --mates2 $fastq_file_2 \
            --geneMap $reference_transcripts_fasta_file \
            --output ${sample_id}_salmon_output/ \
            --threads ${task.cpus}
        """
}


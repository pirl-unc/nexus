#!/usr/bin/env nextflow

// process runLumpyExpressSingleSample {
//
//     label 'lumpy'
//     tag "${sample_id}"
//     debug true
//
//     publishDir(
//         path: "${output_dir}/",
//         mode: 'copy'
//     )
//
//     input:
//         tuple val(sample_id), path(bam_file), path(bam_bai_file)
//         val(python2)
//         val(lumpyexpress)
//         val(lumpyexpress_config_file)
//         val(lumpy_extract_split_reads_script_file)
//         val(samtools)
//         val(svtyper)
//         val(svtyper_sso)
//         val(reference_genome_fasta_file)
//         val(batch_size)
//         val(max_reads)
//         val(output_dir)
//
//     output:
//         tuple val(sample_id), path("${bam_file.baseName}_lumpy.vcf"), path("${bam_file.baseName}_lumpy.gt.vcf"), emit: f
//
//     script:
//         """
//         PYTHONPATH=$python2
//         export PYTHONPATH
//         $samtools view -b -F 1294 $bam_file -o ${bam_file.baseName}_discordants.bam
//         $samtools sort -@ ${task.samtools_cpus} -m ${task.samtools_memory.toGiga()}G ${bam_file.baseName}_discordants.bam -o ${bam_file.baseName}_discordants_sorted.bam
//         $samtools index -b ${bam_file.baseName}_discordants_sorted.bam ${bam_file.baseName}_discordants_sorted.bam.bai
//         $samtools view -h $bam_file | $python2 $lumpy_extract_split_reads_script_file -i stdin | $samtools view -Sb - -o ${bam_file.baseName}_splitters.bam
//         $samtools sort -@ ${task.samtools_cpus} -m ${task.samtools_memory.toGiga()}G ${bam_file.baseName}_splitters.bam -o ${bam_file.baseName}_splitters_sorted.bam
//         $samtools index -b ${bam_file.baseName}_splitters_sorted.bam ${bam_file.baseName}_splitters_sorted.bam.bai
//         $lumpyexpress \
//             -B $bam_file \
//             -S ${bam_file.baseName}_splitters_sorted.bam \
//             -D ${bam_file.baseName}_discordants_sorted.bam \
//             -K $lumpyexpress_config_file \
//             -o ${bam_file.baseName}_lumpy.vcf
//         $svtyper \
//             -B $bam_file \
//             -l ${bam_file.baseName}_svtyper.json
//         $svtyper_sso \
//             --core ${task.cpus} \
//             --batch_size $batch_size \
//             --max_reads $max_reads \
//             -i ${bam_file.baseName}_lumpy.vcf \
//             -B $bam_file \
//             -l ${bam_file.baseName}_svtyper.json \
//             --output_vcf ${bam_file.baseName}_lumpy.gt.vcf
//         """
// }

process runLumpyExpressTumorNormal {

    label 'lumpy'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_lumpy.vcf"), emit: f

    script:
        """
        # Extract the discordant paired-end alignments.
        samtools view -b -F 1294 $tumor_bam_file > ${tumor_bam_file.baseName}_discordants.bam
        samtools view -b -F 1294 $normal_bam_file > ${normal_bam_file.baseName}_discordants.bam
        samtools sort -@ ${task.samtools_cpus} -m ${task.samtools_memory.toGiga()}G ${tumor_bam_file.baseName}_discordants.bam -o ${tumor_bam_file.baseName}_discordants_sorted.bam
        samtools sort -@ ${task.samtools_cpus} -m ${task.samtools_memory.toGiga()}G ${normal_bam_file.baseName}_discordants.bam -o ${normal_bam_file.baseName}_discordants_sorted.bam
        samtools index -b ${tumor_bam_file.baseName}_discordants_sorted.bam ${tumor_bam_file.baseName}_discordants_sorted.bam.bai
        samtools index -b ${normal_bam_file.baseName}_discordants_sorted.bam ${normal_bam_file.baseName}_discordants_sorted.bam.bai

        # Extract the split-read alignments
        samtools view -h $tumor_bam_file | /opt/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - -o ${tumor_bam_file.baseName}_splitters.bam
        samtools view -h $normal_bam_file | /opt/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - -o ${normal_bam_file.baseName}_splitters.bam
        samtools sort -@ ${task.samtools_cpus} -m ${task.samtools_memory.toGiga()}G ${tumor_bam_file.baseName}_splitters.bam -o ${tumor_bam_file.baseName}_splitters_sorted.bam
        samtools sort -@ ${task.samtools_cpus} -m ${task.samtools_memory.toGiga()}G ${normal_bam_file.baseName}_splitters.bam -o ${normal_bam_file.baseName}_splitters_sorted.bam
        samtools index -b ${tumor_bam_file.baseName}_splitters_sorted.bam ${tumor_bam_file.baseName}_splitters_sorted.bam.bai
        samtools index -b ${normal_bam_file.baseName}_splitters_sorted.bam ${normal_bam_file.baseName}_splitters_sorted.bam.bai

        # Run lumpyexpress
        /opt/lumpy-sv/bin/lumpyexpress \
            -B $tumor_bam_file,$normal_bam_file \
            -S ${tumor_bam_file.baseName}_splitters_sorted.bam,${normal_bam_file.baseName}_splitters_sorted.bam \
            -D ${tumor_bam_file.baseName}_discordants_sorted.bam,${normal_bam_file.baseName}_discordants_sorted.bam \
            -o ${sample_id}_lumpy.vcf
        """
}
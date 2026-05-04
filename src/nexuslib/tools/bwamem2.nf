#!/usr/bin/env nextflow

process runBwaMem2Index {

    label 'bwamem2_index'
    debug true

    input:
        path(fasta_file)

    output:
        path("bwamem2_index"), emit: f

    script:
        """
        mkdir -p bwamem2_index
        bwa-mem2 index -p bwamem2_index/$fasta_file $fasta_file
        """
}

process runBwaMem2 {

    label 'bwamem2'
    tag "${sample_id}"
    debug true

    input:
        // fastq_files_1 / fastq_files_2 may each be a single path or a list
        // of paths. When multiple rows in the input TSV share a sample_id,
        // groupTuple in the entry workflow yields parallel R1 / R2 lists
        // here. We cat them in-process (gzip streams concatenate cleanly,
        // and uncompressed FASTQs concatenate trivially) to feed a single
        // bwa-mem2 invocation per sample, producing one merged BAM.
        tuple val(sample_id), path(fastq_files_1, stageAs: 'r1_in/*'), path(fastq_files_2, stageAs: 'r2_in/*')
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_dict_file)
        path(reference_genome_fasta_0123_file)
        path(reference_genome_fasta_amb_file)
        path(reference_genome_fasta_ann_file)
        path(reference_genome_fasta_bwt_file)
        path(reference_genome_fasta_pac_file)
        val(platform_tag)
        val(platform_unit_tag)
        val(library_tag)

    output:
        tuple val(sample_id), path("${sample_id}_bwamem2_sorted.bam"), path("${sample_id}_bwamem2_sorted.bam.bai"), emit: f

    script:
        """
        # Stream the concatenation of all R1s (and all R2s) into bwa-mem2
        # via bash process substitution. No merged FASTQ ever lands on
        # disk — bytes flow cat -> /dev/fd/N -> bwa-mem2 directly.
        #
        # Compression detection is not needed: bwa-mem2 (via kseq.h)
        # detects gzip from magic bytes at the start of each input
        # stream, so plain or gzipped inputs both work transparently.
        # Concatenating gzipped FASTQs is also valid: cat a.fq.gz
        # b.fq.gz produces a multi-member gzip stream that decompresses
        # to the concatenation of a and b.
        #
        # Glob order is alphabetical, so r1_in/* and r2_in/* must use
        # parallel naming conventions for mate pairs to align (the same
        # constraint applied to the previous merged-file approach).
        bwa-mem2 mem -t ${task.bwamem2_threads} \\
            -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:${platform_tag}\\tLB:${library_tag}\\tPU:${platform_unit_tag}" \\
            ${reference_genome_fasta_file} \\
            <(cat r1_in/*) \\
            <(cat r2_in/*) \\
            | samtools view -@ ${task.samtools_view_threads} -bS \\
            | samtools sort -@ ${task.samtools_sort_threads} -m ${task.samtools_memory.toGiga()}G -O bam -o ${sample_id}_bwamem2_sorted.bam
        samtools index -@ ${task.cpus} -b ${sample_id}_bwamem2_sorted.bam ${sample_id}_bwamem2_sorted.bam.bai
        """
}

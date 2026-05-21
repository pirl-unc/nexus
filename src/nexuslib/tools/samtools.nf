#!/usr/bin/env nextflow

process runSamtoolsFaidxFasta {

    label 'samtools_faidx'
    debug true

    input:
        path(fasta_file)

    output:
        path("${fasta_file}"), emit: fasta
        path("${fasta_file}.fai"), emit: fai_file

    script:
        """
        samtools faidx ${fasta_file}
        """
}

process runSamtoolsFaidx {

    label 'samtools_faidx'
    debug true

    input:
        path(fasta_file)

    output:
        path("${gz_name}"), emit: fasta
        path("${gz_name}.fai"), emit: fai_file
        path("${gz_name}.gzi"), emit: gzi_file

    script:
        gz_name = fasta_file.name.endsWith('.gz') ? fasta_file.name : "${fasta_file.name}.gz"
        """
        if [[ "$fasta_file" == *.gz ]]; then
            samtools faidx ${fasta_file}
        else
            bgzip -f $fasta_file
            samtools faidx ${fasta_file}.gz
        fi
        """
}

process runSamtoolsSamToBam {

    label 'samtools_samtobam'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(sam_file)

    output:
        tuple val(sample_id), path("${sam_file.baseName}.bam"), path("${sam_file.baseName}.bam.bai"), emit: f

    script:
        """
        samtools sort -@ ${task.cpus} -m ${task.samtools_memory.toGiga()}G -O bam -o ${sam_file.baseName}.bam $sam_file
        samtools index -@ ${task.cpus} -b ${sam_file.baseName}.bam ${sam_file.baseName}.bam.bai
        """
}

process runSamtoolsSamToBamCustomReference {

    label 'samtools_samtobam_custom_reference'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(sam_file), path(reference_genome_fasta_file), path(reference_genome_fasta_fai_file)

    output:
        tuple val(sample_id), path("${sam_file.baseName}.bam"), path("${sam_file.baseName}.bam.bai"), path(reference_genome_fasta_file), path(reference_genome_fasta_fai_file), emit: f

    script:
        """
        samtools sort -@ ${task.cpus} -m ${task.samtools_memory.toGiga()}G -O bam -o ${sam_file.baseName}.bam $sam_file
        samtools index -@ ${task.cpus} -b ${sam_file.baseName}.bam ${sam_file.baseName}.bam.bai
        """
}

process runSamtoolsSortByName {

    label 'samtools_sort'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_sortedbyname.bam"), emit: f

    script:
        """
        samtools sort -n -@ ${task.cpus} -m ${task.samtools_memory.toGiga()}G -O bam -o ${bam_file.baseName}_sortedbyname.bam $bam_file
        """
}


process runSamtoolsSort {

    label 'samtools_sort'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_sorted.bam"), path("${bam_file.baseName}_sorted.bam.bai"), emit: f

    script:
        """
        samtools sort -@ ${task.cpus} -m ${task.samtools_memory.toGiga()}G -O bam -o ${bam_file.baseName}_sorted.bam $bam_file
        samtools index -@ ${task.cpus} -b ${bam_file.baseName}_sorted.bam ${bam_file.baseName}_sorted.bam.bai
        """
}

process runSamtoolsSortCustomReference {

    label 'samtools_sort_custom_reference'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_sorted.bam"), path("${bam_file.baseName}_sorted.bam.bai"), emit: f

    script:
        """
        samtools sort -@ ${task.cpus} -m ${task.samtools_memory.toGiga()}G -O bam -o ${bam_file.baseName}_sorted.bam $bam_file
        samtools index -@ ${task.cpus} -b ${bam_file.baseName}_sorted.bam ${bam_file.baseName}_sorted.bam.bai
        """
}

process runSamtoolsIndex {

    label 'samtools_index'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file)

    output:
        tuple val(sample_id), path("${bam_file}"), path("${bam_file.baseName}.bam.bai"), emit: f

    script:
        """
        samtools index -@ ${task.cpus} -b ${bam_file} ${bam_file.baseName}.bam.bai
        """
}

process runSamtoolsCalmd {

    label 'samtools_calmd'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_mdtagged.bam"), path("${bam_file.baseName}_mdtagged.bam.bai"), emit: f

    script:
        """
        samtools calmd -@ ${task.cpus} -b $bam_file $reference_genome_fasta_file > ${bam_file.baseName}_mdtagged.bam
        samtools index -@ ${task.cpus} -b ${bam_file.baseName}_mdtagged.bam ${bam_file.baseName}_mdtagged.bam.bai
        """
}

process runSamtoolsCalmdCustomReference {

    label 'samtools_calmd_custom_reference'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(reference_genome_fasta_file), path(reference_genome_fasta_fai_file)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_mdtagged.bam"), path("${bam_file.baseName}_mdtagged.bam.bai"), emit: f

    script:
        """
        ref_file=$reference_genome_fasta_file
        if [[ "$reference_genome_fasta_file" == *.gz ]]; then
            gunzip -c $reference_genome_fasta_file > ref_uncompressed.fa
            samtools faidx ref_uncompressed.fa
            ref_file=ref_uncompressed.fa
        fi
        samtools calmd -@ ${task.cpus} -b $bam_file \$ref_file > ${bam_file.baseName}_mdtagged.bam
        samtools index -@ ${task.cpus} -b ${bam_file.baseName}_mdtagged.bam ${bam_file.baseName}_mdtagged.bam.bai
        """
}

process runSamtoolsMarkdup {

    label 'samtools_markdup'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_markeddup.bam"), path("${bam_file.baseName}_markeddup.bam.bai"), emit: f

    script:
        """
        samtools sort -@ ${task.cpus} -m ${task.samtools_memory.toGiga()}G $bam_file -o ${bam_file.baseName}_position_sorted.bam
        samtools markdup -@ ${task.cpus} ${bam_file.baseName}_position_sorted.bam ${bam_file.baseName}_markeddup.bam
        samtools index -@ ${task.cpus} -b ${bam_file.baseName}_markeddup.bam ${bam_file.baseName}_markeddup.bam.bai
        """
}

process runSamtoolsFixmate {

    label 'samtools_fixmate'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_fixmate.bam"), emit: f

    script:
        """
        samtools sort -n -@ ${task.cpus} -m ${task.samtools_memory.toGiga()}G $bam_file -o ${bam_file.baseName}_name_sorted.bam
        samtools fixmate -m --threads ${task.cpus} ${bam_file.baseName}_name_sorted.bam ${bam_file.baseName}_fixmate.bam
        """
}

process runSamtoolsCoverage {

    label 'samtools_coverage'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(min_mapping_quality)
        val(min_base_quality)
        val(output_dir)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_samtools_coverage.txt"), emit: f

    script:
        """
        samtools coverage \
            --min-MQ ${min_mapping_quality} \
            --min-BQ ${min_base_quality} \
            $bam_file \
            --output ${bam_file.baseName}_samtools_coverage.txt
        """
}

process runSamtoolsFilter {

    label 'samtools_filter'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file)
        val(params_samtools_view)

    output:
        tuple val(sample_id), path("${bam_file.baseName}_samtools_filtered.bam"), emit: f

    script:
        """
        samtools view \
            $params_samtools_view \
            $bam_file \
            --output ${bam_file.baseName}_samtools_filtered.bam
        """
}

process runSamtoolsFastq {

    label 'samtools_fastq'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file)
        val(output_dir)

    output:
        tuple val(sample_id), path("${bam_file.baseName}.fastq.gz"), emit: f

    script:
        """
        samtools fastq --threads ${task.cpus} $bam_file > ${bam_file.baseName}.fastq
        gzip ${bam_file.baseName}.fastq
        """
}

process runSamtoolsMerge {

    label 'samtools_merge'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_files)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_merged.bam"), emit: merged_bam

    script:
        """
        samtools merge -@ ${task.cpus} -o ${sample_id}_merged.bam $bam_files
        """
}


process runSamtoolsMergeRecalibratedBamFiles {

    // Gather half of the scatter-gather ApplyBQSR pattern. Combines N
    // coordinate-sorted, region-disjoint BAM shards into one coordinate-
    // sorted BAM and indexes it. samtools merge re-orders by coordinate
    // using the @SQ header order, so the upstream scatter does not need
    // to emit shards in any particular order.
    //
    //   -c   collapse duplicate @RG header lines across shards
    //   -p   collapse duplicate @PG (program) header lines — otherwise
    //        each shard's ApplyBQSR PG record would appear N times
    //   -f   overwrite any existing output

    label 'samtools_merge'
    tag "${sample_id}"
    debug true

    // bam_basename carries the original (pre-scatter) BAM's basename, so the
    // merged output keeps the same long-form name the unscattered ApplyBQSR
    // would have produced (e.g. "<sample>_bwamem2_sorted_fixmate_markeddup_recalibrated.bam").
    input:
        tuple val(sample_id), val(bam_basename), path(bam_files)

    output:
        tuple val(sample_id), path("${bam_basename}_recalibrated.bam"), path("${bam_basename}_recalibrated.bam.bai"), emit: f

    script:
        """
        samtools merge -@ ${task.cpus} -c -p -f \
            ${bam_basename}_recalibrated.bam \
            ${bam_files}
        samtools index -@ ${task.cpus} -b \
            ${bam_basename}_recalibrated.bam \
            ${bam_basename}_recalibrated.bam.bai
        """
}
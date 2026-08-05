#!/usr/bin/env nextflow

process runIsolaser {

    label 'isolaser'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bam_file), path(bam_bai_file), path(fastq_file), path(gtf_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_dict_file)
        val(params_isolaser)
        val(output_dir)

    output:
        tuple val(sample_id), path("isolaser/"), emit: f

    script:
        """
        mkdir -p isolaser/
        mkdir -p temp/

        # Step 1. isoLASER reads the GTF via pysam.TabixFile, which needs a
        # coordinate-sorted, bgzip-compressed, tabix-indexed GTF. The input GTF
        # may be plain or (b)gzipped, so normalize it to plain text first.
        if [[ "$gtf_file" == *.gz ]]; then
            zcat $gtf_file > input.gtf
        else
            cp $gtf_file input.gtf
        fi

        # `grep "^#" || true`: a header-less GTF (e.g. gencode) makes grep exit 1,
        # which would abort this step under Nextflow's `set -e`.
        (grep "^#" input.gtf || true; grep -v "^#" input.gtf | sort -k1,1 -k4,4n) \
            | bgzip > isolaser/annotation.sorted.gtf.gz
        tabix -p gff isolaser/annotation.sorted.gtf.gz

        # Step 2. Generate transcriptome reference from the indexed GTF
        isolaser_convert_gtf_to_fasta \
            -g isolaser/annotation.sorted.gtf.gz \
            -f $reference_genome_fasta_file \
            -o isolaser/${sample_id}_isolaser.fa

        # Step 3. Align against newly generated transcriptome reference
        minimap2 \
            -t ${task.cpus} \
            -ax splice:hq -uf --MD \
            isolaser/${sample_id}_isolaser.fa \
            $fastq_file > temp/${sample_id}_transcriptome.sam

        # Step 4. Annotate and filter BAM file
        isolaser_annotate \
            -b $bam_file \
            -t temp/${sample_id}_transcriptome.sam \
            -g isolaser/annotation.sorted.gtf.gz \
            -o temp/${sample_id}_transcriptome.bam

        # Step 5. Sort and index the annotated BAM file.
        # isolaser_annotate emits reads in annotation order, so coordinate-sort
        # before indexing (samtools index requires a coordinate-sorted BAM).
        samtools sort -@ ${task.cpus} \
            -o isolaser/${sample_id}_transcriptome.sorted.bam \
            temp/${sample_id}_transcriptome.bam
        samtools index -@ ${task.cpus} -b \
            isolaser/${sample_id}_transcriptome.sorted.bam \
            isolaser/${sample_id}_transcriptome.sorted.bam.bai

        # Step 6. Extract exonic parts from GTF file
        mkdir -p isolaser/transcripts_db/
        isolaser_extract_exon_parts \
            -g isolaser/annotation.sorted.gtf.gz \
            -o isolaser/transcripts_db/

        # Step 7. Run Isolaser
        isolaser \
            -b isolaser/${sample_id}_transcriptome.sorted.bam \
            -o isolaser/${sample_id} \
            -f $reference_genome_fasta_file \
            -t isolaser/transcripts_db/ \
            -s $sample_id \
            $params_isolaser
        """
}

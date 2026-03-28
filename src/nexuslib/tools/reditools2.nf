#!/usr/bin/env nextflow

process prepareGencodeGtfFileForReditools2 {

    label 'reditools2'
    debug true

    input:
        path(gtf_file)

    output:
        path("${gtf_file.baseName}.filtered.sorted.gtf.gz"), emit: gtf_file
        path("${gtf_file.baseName}.filtered.sorted.gtf.gz.tbi"), emit: tbi_file

    script:
        """
        cat $gtf_file | grep -v '^#' | grep 'transcript_id' | sort -k1,1 -k4,4n | bgzip > ${gtf_file.baseName}.filtered.sorted.gtf.gz
        tabix -p gff ${gtf_file.baseName}.filtered.sorted.gtf.gz
        """
}

process runReditools2 {

    label 'reditools2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(rna_bam_file), path(rna_bam_bai_file), path(dna_bam_file), path(dna_bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_genome_fasta_gzi_file)
        path(reference_gtf_file)
        val(params_reditools2)
        val(params_reditools_annotatetable)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_rna_reditools2_final_annotated.tsv"), emit: f

    script:
        """
        # Step 1. Run REDItools2 for the RNA BAM (no header)
        reditools.py \
            -f $rna_bam_file \
            -o ${sample_id}_rna_redtiools2_noheader.tsv \
            -r $reference_genome_fasta_file \
            -H \
            -S \
            $params_reditools2

        # Step 2. Run REDItools2 for the RNA BAM (with header)
        reditools.py \
            -f $rna_bam_file \
            -o ${sample_id}_rna_redtiools2.tsv \
            -r $reference_genome_fasta_file \
            -S \
            $params_reditools2

        # Step 3. Convert the RNA TSV output to BED file
        reditools_table_to_bed.py \
            -i ${sample_id}_rna_redtiools2_noheader.tsv \
            -o ${sample_id}_rna_redtiools2.bed

        # Step 4. Run REDItools2 for the DNA BAM (with header)
        reditools.py \
            -f $dna_bam_file \
            -o ${sample_id}_dna_reditools.tsv \
            -r $reference_genome_fasta_file \
            --dna \
            -B ${sample_id}_rna_redtiools2.bed \
            $params_reditools2

        # Step 5. Identify RNA-specific variants
        annotate_with_DNA.py \
            -r ${sample_id}_rna_redtiools2.tsv \
            -d ${sample_id}_dna_reditools.tsv \
            -R $reference_genome_fasta_fai_file \
            -Z > ${sample_id}_rna_reditools2_final.tsv

        # Step 6. Annotate
        AnnotateTable.py \
            -i ${sample_id}_rna_reditools2_final.tsv \
            -a $reference_gtf_file \
            $params_reditools_annotatetable \
            -o ${sample_id}_rna_reditools2_final_annotated.tsv
        """
}

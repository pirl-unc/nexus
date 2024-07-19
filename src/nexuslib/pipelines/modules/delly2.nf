#!/usr/bin/env nextflow

// process runDelly2SingleSample {
//
//     label 'delly2'
//     tag "${sample_id}"
//     debug true
//
//     publishDir(
//         path: "${output_dir}/",
//         mode: 'symlink'
//     )
//
//     input:
//         tuple val(sample_id), path(bam_file), path(bam_bai_file)
//         val(delly2)
//         val(bcftools)
//         val(reference_genome_fasta_file)
//         val(exclude_regions_file)
//         val(minimum_mapping_quality)
//         val(output_dir)
//
//     output:
//         tuple val(sample_id), path("${bam_file.baseName}_delly2.vcf"), path("${bam_file.baseName}_delly2.bcf"), emit: f
//
//     script:
//         """
//         export OMP_NUM_THREADS=${task.cpus}
//         $delly2 call \
//             --exclude $exclude_regions_file \
//             --outfile ${bam_file.baseName}_delly2.bcf \
//             --genome $reference_genome_fasta_file \
//             --map-qual $minimum_mapping_quality \
//             $bam_file
//         $bcftools view ${bam_file.baseName}_delly2.bcf > ${bam_file.baseName}_delly2.vcf
//         """
// }

process runDelly2TumorNormal {

    label 'delly2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(exclude_tsv_file)
        val(params_delly2call)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_delly2_filtered.vcf"), path("${sample_id}_delly2_filtered.bcf"), emit: f

    script:
        """
        touch delly2_samples.tsv
        tumor_sample_id=\$(samtools view -H "$tumor_bam_file" | grep '^@RG' | awk -F'\t' '{for(i=1;i<=NF;i++) if(\$i ~ /^SM:/) print substr(\$i, 4)}')
        normal_sample_id=\$(samtools view -H "$normal_bam_file" | grep '^@RG' | awk -F'\t' '{for(i=1;i<=NF;i++) if(\$i ~ /^SM:/) print substr(\$i, 4)}')
        printf "%s\t%s\n" "\$tumor_sample_id" "tumor" >> delly2_samples.tsv
        printf "%s\t%s\n" "\$normal_sample_id" "control" >> delly2_samples.tsv
        export OMP_NUM_THREADS=${task.cpus}
        delly call \
            --outfile ${sample_id}_delly2.bcf \
            --genome $reference_genome_fasta_file \
            --exclude $exclude_tsv_file \
            $params_delly2call \
            $tumor_bam_file $normal_bam_file
        delly filter \
            -f somatic \
            -o ${sample_id}_delly2_filtered.bcf \
            -s delly2_samples.tsv \
            ${sample_id}_delly2.bcf
        bcftools view ${sample_id}_delly2_filtered.bcf > ${sample_id}_delly2_filtered.vcf
        """
}
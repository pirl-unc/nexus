#!/usr/bin/env nextflow

process runCircExplorer2Parse {

    label 'circexplorer2'
    tag "${sample_id}"
    debug true

    input:
        tuple val(sample_id), path(input_file)
        val(params_circexplorer2_parse)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_circexplorer2_back_spliced_junction.bed"), emit: f

    script:
        """
        # Check if first line is a header (non-numeric in coordinate column)
        if head -1 $input_file | awk '{exit (\$2+0 == \$2) ? 1 : 0}'; then
            tail -n +2 $input_file > input_no_header.junction
        else
            ln -s $input_file input_no_header.junction
        fi

        CIRCexplorer2 parse \
            -b ${sample_id}_circexplorer2_back_spliced_junction.bed \
            $params_circexplorer2_parse \
            input_no_header.junction
        """
}

process runCircExplorer2Annotate {

    label 'circexplorer2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(back_spliced_junction_bed_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(circexplorer2_gene_annotation_txt_file)
        val(params_circexplorer2_annotate)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_circexplorer2_back_spliced_junction_annotated.bed"), emit: f

    script:
        """
        CIRCexplorer2 annotate \
            $params_circexplorer2_annotate \
            -r $circexplorer2_gene_annotation_txt_file \
            -g $reference_genome_fasta_file \
            -b $back_spliced_junction_bed_file \
            -o ${sample_id}_circexplorer2_back_spliced_junction_annotated.bed
        """
}

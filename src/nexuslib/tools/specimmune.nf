#!/usr/bin/env nextflow

process runSpecImmune {

    label 'specimmune'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy',
        // Strip the "${sample_id}_specimmune_outputs/${sample_id}/" prefix so
        // each child entry lands directly under ${output_dir}/${sample_id}/.
        // Without this, publishDir preserves the relative path from the work
        // dir and you end up with
        //   ${output_dir}/${sample_id}/${sample_id}_specimmune_outputs/${sample_id}/<file>
        saveAs: { f ->
            def prefix = "${sample_id}_specimmune_outputs/${sample_id}/"
            f.startsWith(prefix) ? f.substring(prefix.length()) : f
        }
    )

    input:
        tuple val(sample_id), path(fastq_file)
        val(params_specimmune)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_specimmune_outputs/${sample_id}/*"), emit: f

    script:
        """
        mkdir -p ${sample_id}_specimmune_outputs/
        python /opt/SpecImmune-0.0.3/scripts/main.py \
            -r $fastq_file \
            -n $sample_id \
            -o ${sample_id}_specimmune_outputs/ \
            --db /opt/SpecImmune-0.0.3/db/ \
            -j ${task.cpus} \
            $params_specimmune
        """
}
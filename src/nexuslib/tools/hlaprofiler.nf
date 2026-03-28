#!/usr/bin/env nextflow

process runHLAProfilerPredict {

    label 'hlaprofiler'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(params_hlaprofiler)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_hlaprofiler_output/*"), emit: f

    script:
        """
        mkdir ${sample_id}_hlaprofiler_output/
        HLAProfiler.pl predict \
          -fastq1 $fastq_file_1 \
          -fastq2 $fastq_file_2 \
          -threads ${task.cpus} \
          -output_dir ${sample_id}_hlaprofiler_output/ \
          -kraken_path /opt/kraken/kraken-ea-0.10.5ea.3-3 \
          -database_dir /opt/HLAProfiler/HLAProfiler-1.0.0-db_only/ \
          -database_name hla_database \
          -reference /opt/HLAProfiler/HLAProfiler-1.0.0-db_only/hla_database/data/reference/hla.ref.merged.fa \
          -l ${sample_id}_hlaprofiler_output/${sample_id}_HLAProfiler.log \
          $params_hlaprofiler
        """
}
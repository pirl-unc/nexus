#!/usr/bin/env nextflow

process runHLAProfilerPredict {

    label 'hlaprofiler'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fastq_file_1), path(fastq_file_2)
        val(params_hlaprofiler)
        val(output_dir)

    output:
        tuple val(sample_id), path("hlaprofiler/"), emit: f

    script:
        """
        mkdir hlaprofiler/
        HLAProfiler.pl predict \
          -fastq1 $fastq_file_1 \
          -fastq2 $fastq_file_2 \
          -threads ${task.cpus} \
          -output_dir hlaprofiler/ \
          -kraken_path /opt/kraken/kraken-ea-0.10.5ea.3-3 \
          -database_dir /opt/HLAProfiler/HLAProfiler-1.0.0-db_only/ \
          -database_name hla_database \
          -reference /opt/HLAProfiler/HLAProfiler-1.0.0-db_only/hla_database/data/reference/hla.ref.merged.fa \
          -l hlaprofiler/${sample_id}_HLAProfiler.log \
          $params_hlaprofiler
        """
}

#!/usr/bin/env nextflow

process runBeers2 {

    label 'beers2'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/beers2/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(beers_config_file)
        val(output_dir)

    output:
        // Emit the three top-level BEERS2 output directories:
        //   library_prep_pipeline/ — per-sample lib-prep packets + logs
        //   sequence_pipeline/      — per-sample flowcell cluster packets + logs
        //   results/                — final FASTQ/SAM/BAM files
        tuple val(sample_id),
              path("library_prep_pipeline/"),
              path("sequence_pipeline/"),
              path("results/"),
              emit: f

    script:
        """
        beers --configfile $beers_config_file --cores ${task.cpus}

        # Gzip every FASTQ in the BEERS2 results directory. xargs -P
        # parallelizes one gzip per file across ${task.cpus}.
        find results/ -type f -name "*.fastq" -print0 \
            | xargs -0 -r -n 1 -P ${task.cpus} gzip
        """
}

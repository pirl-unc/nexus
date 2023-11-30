## long_read_alignment_minimap2.nf

Aligns long DNA or RNA reads (Pacific Biosciences or Oxford Nanopore) to a reference genome using [minimap2](https://github.com/lh3/minimap2).

### Inputs / Outputs

| I/O    | Description                                                      |
|:-------|:-----------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                                 | 
| Output | MD-tagged and sorted `bam` and `bam.bai` files for each sample.  |

### Dependencies

* `minimap2`
* `samtools`

### Usage

```
workflow:
    1. Align reads (fastq.gz files) to a reference genome using minimap2.
    2. Convert sam files to bam files.
    3. Generate MD tags.
    4. Sort MD-tagged bam files.

usage: nexus run --nf-workflow long_read_alignment_minimap2.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
    --output_dir                    :   Directory to which output files will be copied.
    --reference_genome_fasta_file   :   Reference genome FASTA file.
    --minimap2                      :   minimap2 path.
    --minimap2_params               :   minimap2 parameters (e.g. '"--ax map-hifi --cs --eqx --Y --L "').
                                        Note that the parameters need to be wrapped in quotes 
                                        and a space at the end of the string is necessary.
    --samtools                      :   samtools path.

optional arguments:
    --platform_tag                  :   Platform tag (default: 'unknown').
    --platform_unit_tag             :   Platform unit tag (default: 'unknown').
    --library_tag                   :   Library tag (default: 'unknown').
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/configs)

`--sample_tsv_file`

| Header     | Description                   |
| ---------- |-------------------------------|
| sample_id  | Sample ID.                    |
| fastq_file | Full path to `fastq.gz` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on the LBG.

`--minimap2_params`
* Refer to the [minimap2 documentation](https://github.com/lh3/minimap2).
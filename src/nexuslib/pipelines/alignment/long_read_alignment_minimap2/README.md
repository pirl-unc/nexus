## long_read_alignment_minimap2.nf

Aligns long DNA or RNA reads (Pacific Biosciences or Oxford Nanopore) to a reference genome using [minimap2](https://github.com/lh3/minimap2).

### Inputs / Outputs

| Input(s)        | Output(s)                                       |
| --------------- | ----------------------------------------------- |
| `FASTQ.GZ` file | Sorted and MD-tagged `BAM` and `BAM.BAI` files  |

### Dependencies

* `minimap2`
* `samtools`

### Usage

```shell
nextflow run long_read_alignment_minimap2.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
    --output_dir                    :   Directory to which output files will be copied.
    --reference_genome_fasta_file   :   Reference genome FASTA file.
    --minimap2                      :   minimap2 path.
    --minimap2_params               :   minimap2 parameters (e.g. '--ax map-hifi --cs --eqx --Y --L').
    --samtools                      :   samtools path.

optional arguments:
    --platform_tag                  :   Platform tag (default: 'unknown').
    --platform_unit_tag             :   Platform unit tag (default: 'unknown').
    --library_tag                   :   Library tag (default: 'unknown').
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`--sample_tsv_file`

| Header     | Description                  |
| ---------- | ---------------------------- |
| sample_id  | Sample ID                    |
| fastq_file | Full path to `FASTQ.GZ` file |

`--minimap2_params`

Refer to the [minimap2 documentation](https://github.com/lh3/minimap2).
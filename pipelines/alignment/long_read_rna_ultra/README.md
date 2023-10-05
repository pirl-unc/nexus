## long_read_rna_alignment_ultra.nf

Aligns long RNA reads (Pacific Biosciences or Oxford Nanopore) to a reference genome using [uLTRA](https://github.com/ksahlin/ultra).

### Inputs / Outputs

| Input(s)        | Output(s)                                       |
| --------------- | ----------------------------------------------- |
| `FASTQ.GZ` file | Sorted and MD-tagged `BAM` and `BAM.BAI` files  |

### Dependencies

* `uLTRA`
* `samtools`

### Usage

```shell
nextflow run long_read_rna_alignment_ultra.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
    --output_dir                    :   Directory to which output files will be copied.
    --reference_genome_fasta_file   :   Reference genome FASTA file.
    --ultra                         :   uLTRA path.
    --ultra_index                   :   uLTRA index path.
    --ultra_params                  :   uLTRA parameters (e.g. '--isoseq').
    --samtools                      :   samtools path.

optional arguments:
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`--sample_tsv_file`

| Header     | Description                |
| ---------- | -------------------------- |
| sample_id  | Sample ID                  |
| fastq_file | Full path to `FASTQ.GZ` file |

`--ultra_index`

An [uLTRA index](https://github.com/ksahlin/ultra) needs to be generated prior to running this workflow script.

`--ultra_params`

Refer to the [uLTRA documentation](https://github.com/ksahlin/ultra).
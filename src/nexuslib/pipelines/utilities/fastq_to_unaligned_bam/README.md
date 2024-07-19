## fastq_to_unaligned_bam.nf

Convert FASTQ files to unaligned BAM files using [Picard](https://broadinstitute.github.io/picard/)

### Inputs / Outputs

| I/O    | Description                           |
|:-------|:--------------------------------------|
| Input  | `fastq.gz` file for each sample.      | 
| Output | Unaligned `bam` file for each sample. |

### Dependencies

* `Picard`
* `samtools`

### Example

```
nexus run --nf-workflow fastq_to_unaligned_bam.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --read_type 'single-end' OR 'paired-end'
```

### Usage

```
workflow:
    1. Run Picard 'FastqToSam'.
    2. Run samtools to convert SAM files to BAM files.

usage: nexus run --nf-workflow fastq_to_unaligned_bam.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --read_type                         :   read type (default: 'single-end'). Either 'single-end' or 'paired-end'.
    --params_picard                     :   Picard 'FastqToSam' parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`

| Header       | Description                                                            |
|--------------|------------------------------------------------------------------------|
| sample_id    | Sample ID.                                                             |
| fastq_file_1 | Full path to `fastq.gz` file.                                          |
| fastq_file_2 | Full path to `fastq.gz` file. (if single-end, leave this column empty) |

`--read_type`
* Either `single-end` or `paired-end`

`--params_picard`
* Refer to the [Picard documentation](https://broadinstitute.github.io/picard/).
* The following parameters for `Picard` are already included in `nexus` module for `picard.jar FastqToSam` and should not be specified:
  * `F1`
  * `F2`
  * `SM`
  * `O`

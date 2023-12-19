## fastq_to_unaligned_bam.nf

Converts FASTQ files to unaligned BAM files using [picard](https://broadinstitute.github.io/picard/)

### Inputs / Outputs

| I/O    | Description                           |
|:-------|:--------------------------------------|
| Input  | `fastq` file for each sample.         | 
| Output | Unaligned `bam` file for each sample. |

### Dependencies

* `picard`
* `samtools`

### Usage

```
workflow:
    1. Run picard 'FastqToSam'.
    2. Run samtools to convert SAM files to BAM files.

usage: nexus run --nf-workflow fastq_to_unaligned_bam.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --read_type                         :   read type (default: 'single-end'). Either 'single-end' or 'paired-end'.
    --samtools                          :   samtools path (default: samtools).
    --picard                            :   picard path (default:
                                            /datastore/lbcfs/collaborations/pirl/share/apps/picard/v2.27.5/picard.jar).
    --picard_params                     :   picard 'FastqToSam' parameters (default: '" "').
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/configs)

`--sample_tsv_file`

| Header       | Description                                                            |
|--------------|------------------------------------------------------------------------|
| sample_id    | Sample ID.                                                             |
| fastq_file_1 | Full path to `fastq.gz` file.                                          |
| fastq_file_2 | Full path to `fastq.gz` file. (if single-end, leave this column empty) |

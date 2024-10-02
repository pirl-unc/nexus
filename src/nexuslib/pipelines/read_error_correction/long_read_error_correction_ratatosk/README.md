## long_read_error_correction_ratatosk.nf

Performs error-correction on long-read sequencing reads using [Ratatosk](https://github.com/DecodeGenetics/Ratatosk).

### Inputs / Outputs

| I/O    | Description                                                                           |
|:-------|:--------------------------------------------------------------------------------------|
| Input  | Long-read `fastq.gz` and short-read (paired-end) `fastq.gz` files for each sample.    | 
| Output | Error-corrected `fastq.gz` file for each sample.                                     |

### Dependencies

* `Ratatosk`

### Example

```
nexus run --nf-workflow long_read_error_correction_ratatosk.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR
```

### Usage

```
workflow:
    1. Correct reads using Ratatosk (1st pass).
    2. Correct reads using Ratatosk (2nd pass).

usage: nexus run --nf-workflow long_read_error_correction_ratatosk.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'long_read_fastq_file', 'short_read_r1_fastq_file', 'short_read_r2_fastq_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --params_ratatosk_first_correct     :   Ratatosk first-pass 'correct' parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
    --params_ratatosk_second_correct    :   Ratatosk second-pass 'correct' parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.

optional arguments:
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header                   | Description                                 |
|--------------------------|---------------------------------------------|
| sample_id                | Sample ID.                                  |
| long_read_fastq_file     | Full path to long-read `fastq.gz` file.     |
| short_read_r1_fastq_file | Full path to short-read R1 `fastq.gz` file. |
| short_read_r2_fastq_file | Full path to short-read R2 `fastq.gz` file. |

`--params_ratatosk_first_correct`
* Refer to the [Ratatosk documentation](https://github.com/DecodeGenetics/Ratatosk).
* The following parameters for `Ratatosk` are already included in `nexus` module for `Ratatosk correct -1` and should not be specified:
  * `-v`
  * `-c`
  * `-l`
  * `-s`
  * `-o`

`--params_ratatosk_second_correct`
* Refer to the [Ratatosk documentation](https://github.com/DecodeGenetics/Ratatosk).
* The following parameters for `Ratatosk` are already included in `nexus` module for `Ratatosk correct -2` and should not be specified:
  * `-v`
  * `-G`
  * `-c`
  * `-l`
  * `-s`
  * `-L`
  * `-o`

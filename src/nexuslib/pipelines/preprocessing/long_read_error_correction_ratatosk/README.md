## long_read_error_correction_ratatosk.nf

Performs error-correction on long-read sequencing data using [ratatosk](https://github.com/DecodeGenetics/Ratatosk).

### Inputs / Outputs

| I/O    | Description                                                                           |
|:-------|:--------------------------------------------------------------------------------------|
| Input  | Long-read `fastq.gz` and short-read (paired-end) `fastq.gz` files for each sample.    | 
| Output | Error-corrected `fastq.gz` file for each sample.                                     |

### Dependencies

* `ratatosk`

### Usage

```
 workflow:
    1. Index using ratatosk (1st pass).
    2. Correct using ratatosk (1st pass).
    3. Index using ratatosk (2nd pass).
    4. Correct using ratatosk (2nd pass).

usage: nexus run --nf-workflow long_read_error_correction_ratatosk.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'long_read_fastq_file', 'short_read_r1_fastq_file', 'short_read_r2_fastq_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --ratatosk                          :   ratatosk path.
    --ratatosk_first_index_params       :   ratatosk first-pass 'index' parameters (default: '" "').
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
    --ratatosk_first_correct_params     :   ratatosk first-pass 'correct' parameters (default: '" "').
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
    --ratatosk_second_index_params      :   ratatosk second-pass 'index' parameters (default: '" "').
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
    --ratatosk_second_correct_params    :   ratatosk second-pass 'correct' parameters (default: '" "').
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.

optional arguments:
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/configs)

`--sample_tsv_file`

| Header     | Description                             |
| ---------- |-----------------------------------------|
| sample_id  | Sample ID.                              |
| long_read_fastq_file | Full path to long-read `fastq.gz` file. |
| short_read_r1_fastq_file | Full path to short-read R1 `fastq.gz` file. |
| short_read_r2_fastq_file | Full path to short-read R2 `fastq.gz` file. |

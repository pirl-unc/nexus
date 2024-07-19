## paired-end_read_rna_hla_typing_seq2hla.nf

Performs HLA typing using paired-end (Illumina) RNA reads using [seq2HLA](https://github.com/TRON-Bioinformatics/seq2HLA).

### Inputs / Outputs

| I/O    | Description                                      |
|:-------|:-------------------------------------------------|
| Input  | `R1 and R2 fastq.gz` files for each sample.      | 
| Output | seq2HLA output files for each sample.            |

### Dependencies

* `seq2HLA`

### Example

```
nexus run --nf-workflow paired-end_read_rna_hla_typing_seq2hla.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR
```

### Usage

```
workflow:
    1. Profile HLA alleles using paired-end RNA sequencing FASTQ files using Seq2HLA.

usage: nexus run --nf-workflow paired-end_read_rna_hla_typing_seq2hla.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --params_seq2hla                :   Seq2HLA parameters (default: '""').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`

| Header       | Description                     |
| ------------ |---------------------------------|
| sample_id    | Sample ID.                      |
| fastq_file_1 | Full path to R1 `fastq.gz` file |
| fastq_file_2 | Full path to R2 `fastq.gz` file |

`params_seq2hla`
* Refer to the [seq2HLA documentation](https://github.com/TRON-Bioinformatics/seq2HLA).
* The following parameters for `seq2HLA` are already included in `nexus` module for `seq2HLA` and should not be specified:
  * `-1`
  * `-2`
  * `-r`
  * `-p`

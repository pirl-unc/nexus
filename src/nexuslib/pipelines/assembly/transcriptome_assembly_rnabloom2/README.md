## transcriptome_assembly_rnabloom2.nf

Assemble transcript sequences using long RNA reads using [RNA-Bloom2](https://github.com/bcgsc/RNA-Bloom).

### Inputs / Outputs

| I/O    | Description                                                   |
|:-------|:--------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                             | 
| Output | RNA-Bloom2 output files including FASTA file for each sample. |

### Dependencies

* `RNA-Bloom2`

### Example

```
nexus run --nf-workflow transcriptome_assembly_rnabloom2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --params_rnabloom2 '"--qual 20 --qual-avg 20 --mincov 3 -ntcard -savebf"'
```

### Usage

```
workflow:
    1. Assemble transcripts using RNA-Bloom2.

usage: nexus run --nf-workflow transcriptome_assembly_rnabloom2.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --params_rnabloom2                  :   RNA-Bloom2 parameters (default: '"--qual 20 --qual-avg 20 --mincov 3 -ntcard -savebf"').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                     |
|--------------|---------------------------------|
| sample_id    | Sample ID.                      |
| fastq_file   | Full path to `fastq.gz` file    |

`--params_rnabloom2`
* Refer to the [RNA-Bloom2 documentation](https://github.com/bcgsc/RNA-Bloom).
* The following parameters for `RNA-Bloom2` are already included in `nexus` module for `RNA-Bloom2` and should not be specified:
  * `-long`
  * `--threads`
  * `--outdir`

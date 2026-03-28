## quantification_kallisto-pe.nf

Quantify RNA abundance in paired-end RNA reads using [Kallisto](https://github.com/pachterlab/kallisto).

### Inputs / Outputs

| I/O    | Description                                                           |
|:-------|:----------------------------------------------------------------------|
| Input  | `R1 and R2 fastq.gz` files for each sample.                           | 
| Output | Kallisto output files including `abundance.tsv` file for each sample. |

### Dependencies

* `Kallisto`

### Example

```
nexus run --nf-workflow quantification_kallisto-pe.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --reference_transcripts_fasta_file REFERENCE_TRANSCRIPTS_FASTA_FILE \
    --output_dir OUTPUT_DIR \
    --params_kallisto_quant '""'
```

### Usage

```
workflow:
    1. Quantify RNA in paired-end read FASTQ files using Kallisto.

usage: nexus run --nf-workflow quantification_kallisto-pe.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --reference_transcripts_fasta_file  :   Reference transcripts FASTA file.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --params_kallisto_quant             :   Kallisto quant parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.

```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                     |
| ------------ |---------------------------------|
| sample_id    | Sample ID.                      |
| fastq_file_1 | Full path to R1 `fastq.gz` file |
| fastq_file_2 | Full path to R2 `fastq.gz` file |

`--reference_transcripts_fasta_file`
* Reference genome FASTA (`.fasta` or `fasta.gz`) file.

* `--params_kallisto_quant`
* Refer to the [Kallisto documentation](https://pachterlab.github.io/kallisto/manual.html).
* The following parameters for `Kallisto` are already included in `nexus` module for `kallisto quant` and should not be specified:
  * `--index`
  * `--output-dir`
  * `--threads`

## quantification_kallisto-lr.nf

Quantify RNA abundance in long RNA reads using [Kallisto](https://github.com/pachterlab/kallisto).

### Inputs / Outputs

| I/O    | Description                                                           |
|:-------|:----------------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                                      | 
| Output | Kallisto output files including `abundance.tsv` file for each sample. |

### Dependencies

* `Kallisto`

### Example

```
nexus run --nf-workflow quantification_kallisto-lr.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --reference_transcripts_fasta_file REFERENCE_TRANSCRIPTS_FASTA_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --output_dir OUTPUT_DIR \
    --params_kallisto_index '"-k 63"' \
    --params_kallisto_bus '"-x bulk --threshold 0.8"' \
    --params_bustools_sort '""' \
    --params_bustools_count '"--cm -m"' \
    --params_kallisto_quanttcc '"-P PacBio --matrix-to-files"'
```

### Usage

```
workflow:
    1. Quantify RNA in long-read FASTQ files using Kallisto (quant-tcc).

usage: nexus run --nf-workflow quantification_kallisto-lr.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file'.
    --reference_transcripts_fasta_file  :   Reference transcripts FASTA file.
    --reference_genes_gtf_file          :   Reference genes GTF file.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --params_kallisto_index             :   Kallisto index parameters (default: '"-k 63"').
    --params_kallisto_bus               :   Kallisto bus parameters (default: '"-x bulk --threshold 0.8"').
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
    --params_bustools_sort              :   Bustools sort parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
    --params_bustools_count             :   Bustools count parameters (default: '"--cm -m"').
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
    --params_kallisto_quanttcc          :   Kallisto quant-tcc parameters (default: '"-P PacBio --matrix-to-files"').
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                   |
|--------------|-------------------------------|
| sample_id    | Sample ID.                    |
| fastq_file   | Full path to `fastq.gz` file. |

`--reference_transcripts_fasta_file`
* Reference genome FASTA (`.fasta` or `fasta.gz`) file.

`--reference_genes_gtf_file`
* Reference genes GTF (`.gtf` or `.gtf.gz`) file.

`--params_kallisto_index`
* Refer to the [Kallisto documentation](https://pachterlab.github.io/kallisto/manual.html).
* The following parameters for `Kallisto` are already included in `nexus` module for `kallisto index` and should not be specified:
  * `--index`
  * `--threads`

`--params_kallisto_bus`
* Refer to the [Kallisto documentation](https://pachterlab.github.io/kallisto/manual.html).
* The following parameters for `Kallisto` are already included in `nexus` module for `kallisto bus` and should not be specified:
  * `-t`
  * `--long`
  * `-i`
  * `-o`

`--params_bustools_sort`
* Refer to the [Kallisto documentation](https://pachterlab.github.io/kallisto/manual.html).
* The following parameters for `bustools` are already included in `nexus` module for `bustools sort` and should not be specified:
  * `-t`
  * `-o`

`--params_bustools_count`
* Refer to the [Kallisto documentation](https://pachterlab.github.io/kallisto/manual.html).
* The following parameters for `bustools` are already included in `nexus` module for `bustools count` and should not be specified:
  * `-t`
  * `-e`
  * `-o`
  * `-g`

`--params_kallisto_quanttcc`
* Refer to the [Kallisto documentation](https://pachterlab.github.io/kallisto/manual.html).
* The following parameters for `Kallisto` are already included in `nexus` module for `kallisto quant-tcc` and should not be specified:
  * `-t`
  * `--long`
  * `-f`
  * `-i`
  * `-e`
  * `-o`
  * `--gtf`

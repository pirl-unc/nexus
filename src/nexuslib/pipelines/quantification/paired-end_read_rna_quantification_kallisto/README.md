## paired-end_read_rna_quantification_kallisto.nf

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
nexus run --nf-workflow paired-end_read_rna_quantification_kallisto.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --kallisto_index_file KALLISTO_INDEX_FILE \
    --output_dir OUTPUT_DIR \
    --params_kallisto_quant PARAMS_KALLISTO_QUANT
```

### Usage

```
workflow:
    1. Quantify RNA in paired-end read FASTQ files using Kallisto.

usage: nexus run --nf-workflow paired-end_read_rna_quantification_kallisto.nf [required] [optional] [--help]

required arguments:
    -c                          :   Nextflow .config file.
    -w                          :   Nextflow work directory path.
    --samples_tsv_file          :   TSV file with the following columns:
                                    'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --kallisto_index_file       :   Kallisto index file.
    --output_dir                :   Directory to which output files will be copied.

optional arguments:
    --params_kallisto_quant     :   Kallisto quant parameters (default: '""').
                                    Note that the parameters need to be wrapped in quotes.
    --delete_work_dir           :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--kallisto_index_file`
* Kallisto index files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/kallisto/` on LBG.

`--params_kallisto_quant`
* Refer to the [Kallisto documentation](https://pachterlab.github.io/kallisto/manual.html).
* The following parameters for `Kallisto` are already included in `nexus` module for `kallisto quant` and should not be specified:
  * `--index`
  * `--output-dir`
  * `--threads`


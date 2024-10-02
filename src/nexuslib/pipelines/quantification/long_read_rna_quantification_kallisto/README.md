## long_read_rna_quantification_kallisto.nf

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
nexus run --nf-workflow long_read_rna_quantification_kallisto.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --kallisto_index_file KALLISTO_INDEX_FILE \
    --output_dir OUTPUT_DIR \
    --params_kallisto_quant_fragment_length PARAMS_KALLISTO_QUANT_FRAGMENTATION_LENGTH \
    --params_kallisto_quant_sd PARAMS_KALLISTO_QUANT_SD \
    --params_kallisto_quant PARAMS_KALLISTO_QUANT
```

### Usage

```
workflow:
    1. Quantify RNA in long-read FASTQ files using Kallisto.

usage: nexus run --nf-workflow long_read_rna_quantification_kallisto.nf [required] [optional] [--help]

required arguments:
    -c                                      :   Nextflow .config file.
    -w                                      :   Nextflow work directory path.
    --samples_tsv_file                      :   TSV file with the following columns:
                                                'sample_id', 'fastq_file'.
    --kallisto_index_file                   :   Kallisto index file.
    --params_kallisto_quant_fragment_length :   Kallisto quant '--fragment-length' parameter value (required for single-end reads).
    --params_kallisto_quant_sd              :   Kallisto quant '--sd' parameter value (required for single-end reads).
    --output_dir                            :   Directory to which output files will be copied.

optional arguments:
    --params_kallisto_quant                 :   Kallisto quant parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                       :   Delete work directory (default: false).
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

`--kallisto_index_file`
* Kallisto index files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/kallisto/` on LBG.

`--params_kallisto_quant`
* Refer to the [Kallisto documentation](https://pachterlab.github.io/kallisto/manual.html).
* The following parameters for `Kallisto` are already included in `nexus` module for `kallisto quant` and should not be specified:
  * `--index`
  * `--output-dir`
  * `--single`
  * `--threads`
  * `--params_kallisto_quant_fragment_length`
  * `--params_kallisto_quant_sd`


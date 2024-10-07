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
    --gtf_file GTF_FILE \
    --t2g_file T2G_FILE \
    --output_dir OUTPUT_DIR \
    --params_kallisto_bus PARAMS_KALLISTO_BUS \
    --params_bustools_sort PARAMS_BUSTOOLS_SORT \
    --params_bustools_count PARAMS_BUSTOOLS_COUNT \
    --params_kallisto_quanttcc PARAMS_KALLISTO_QUANTTCC
```

### Usage

```
workflow:
    1. Quantify RNA in long-read FASTQ files using Kallisto (quant-tcc).

usage: nexus run --nf-workflow long_read_rna_quantification_kallisto.nf [required] [optional] [--help]

required arguments:
    -c                                      :   Nextflow .config file.
    -w                                      :   Nextflow work directory path.
    --samples_tsv_file                      :   TSV file with the following columns:
                                                'sample_id', 'fastq_file'.
    --kallisto_index_file                   :   Kallisto index file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/kallisto/v0.51.1/gencode_v41/kallisto_gencode_v41_index_k63.idx).
    --gtf_file                              :   GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf).
    --t2g_file                              :   T2G file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/kallisto/v0.51.1/gencode_v41/gencode_v41.t2g).
    --params_kallisto_bus                   :   Kallisto bus parameters (default: '"-x bulk --threshold 0.8"').
    --params_bustools_sort                  :   Bustools sort parameters (default: '""').
    --params_bustools_count                 :   Bustools count parameters (default: '"--cm -m"').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
    --params_kallisto_quanttcc              :   Kallisto quant-tcc parameters (default: '"-P PacBio"').
    --output_dir                            :   Directory to which output files will be copied.

optional arguments:
    --delete_work_dir                       :   Delete work directory (default: false).Z
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

`gtf_file`
* GTF files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

`t2g_file`
* T2G files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/kallisto/` on LBG.

`params_kallisto_bus`
* Refer to the [Kallisto documentation](https://pachterlab.github.io/kallisto/manual.html).
* The following parameters for `Kallisto` are already included in `nexus` module for `kallisto bus` and should not be specified:
  * `-t`
  * `--long`
  * `-i`
  * `-o`

`params_bustools_sort`
* Refer to the [Kallisto documentation](https://pachterlab.github.io/kallisto/manual.html).
* The following parameters for `bustools` are already included in `nexus` module for `bustools sort` and should not be specified:
  * `-t`
  * `-o`

`params_bustools_count`
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

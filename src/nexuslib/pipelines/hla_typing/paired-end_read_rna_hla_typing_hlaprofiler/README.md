## paired-end_read_rna_hla_typing_hlaprofiler.nf

Performs HLA typing using paired-end (Illumina) RNA reads using [HLAProfiler](https://github.com/ExpressionAnalysis/HLAProfiler).

### Inputs / Outputs

| I/O    | Description                                                |
|:-------|:-----------------------------------------------------------|
| Input  | `R1 and R2 fastq.gz` files for each sample.       | 
| Output | HLAProfiler outputs including `txt` files for each sample. |

### Dependencies

* `HLAProfiler`

### Example

```
nexus run --nf-workflow paired-end_read_rna_hla_typing_hlaprofiler.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --params_hlaprofiler '"-allele_refinement all -if"'
```

### Usage

```
workflow:
    1. Profile HLA alleles using paired-end RNA sequencing FASTQ files using HLAProfiler.

usage: nexus run --nf-workflow paired-end_read_rna_hla_typing_hlaprofiler.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --params_hlaprofiler            :   HLAProfiler 'predict' parameters (default: '"-allele_refinement all -if"').
                                        Note that the parameters need to be wrapped in quotes.
    --delete_work_dir               :   Delete work directory (default: false).
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

`params_hlaprofiler`
* Refer to the [HLAProfiler documentation](https://github.com/ExpressionAnalysis/HLAProfiler).
* The following parameters for `HLAProfiler` are already included in `nexus` module for `HLAProfiler` and should not be specified:
  * `-fastq1`
  * `-fastq2`
  * `-threads`
  * `-output_dir`
  * `-kraken_path`
  * `-database_dir`
  * `-database_name`
  * `-reference`
  * `-l`

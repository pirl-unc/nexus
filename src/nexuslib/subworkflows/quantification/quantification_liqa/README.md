## quantification_liqa.nf

Quantify RNA abundance in long RNA reads using [LIQA](https://github.com/WGLab/LIQA).

### Inputs / Outputs

| I/O    | Description                                |
|:-------|:-------------------------------------------|
| Input  | `bam` and `bam.bai` files for each sample. | 
| Output | LIQA output files                          |

### Dependencies

* `LIQA`

### Example

```
nexus run --nf-workflow quantification_liqa.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --output_dir OUTPUT_DIR \
    --params_liqa_quantify '"-max_distance 20 -f_weight 1"'
```

### Usage

```
workflow:
    1. Quantify RNA in long-read FASTQ files using LIQA.

usage: nexus run --nf-workflow quantification_liqa.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'bam_file', 'bam_bai_file'.
    --reference_genes_gtf_file      :   Reference genes GTF file.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --params_liqa_quantify          :   LIQA quantify parameters (default: '"-max_distance 20 -f_weight 1"').
                                        Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                  |
|--------------|------------------------------|
| sample_id    | Sample ID.                   |
| bam_file     | Full path to `bam` file.     |
| bam_bai_file | Full path to `bam.bai` file. |

`--reference_genes_gtf_file`
* Reference genes GTF (`.gtf` or `.gtf.gz`) file.

`--params_liqa_quantify`
* Refer to the [LIQA documentation](https://github.com/WGLab/LIQA/blob/master/doc/Usage.md).
* The following parameters for `LIQA` are already included in `nexus` module for `liqa -task quantify` and should not be specified:
  * `-refgene`
  * `-bam`
  * `-out`

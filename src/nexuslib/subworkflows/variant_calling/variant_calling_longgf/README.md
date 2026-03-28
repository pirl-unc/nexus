## variant_calling_longgf.nf

Identifies fusion genes in long-read RNA BAM files using [LongGF](https://github.com/WGLab/LongGF).

### Inputs / Outputs

| I/O    | Description                               |
|:-------|:------------------------------------------|
| Input  | `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.               |

### Dependencies

* `LongGF`

### Example

```
nexus run --nf-workflow variant_calling_longgf.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --params_longgf '"100 30 100 1 0 3"'
```

### Usage

```
workflow:
    1. Run Longgf.

usage: nexus run --nf-workflow variant_calling_longgf.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns: 'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                    :   Directory to which output files will be copied.
    --reference_genes_gtf_file      :   Reference genes annotation GTF file.

optional arguments:
    --params_longgf                 :   Longgf parameters (default: '"100 30 100 1 0 3"').
                                        <min-overlap-len> <bin_size> <min-map-len> [pseudogene:0(default)/1/other(no filter)] [Secondary_alignment:0(default)] [min_sup_read:2(default)]
                                        Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header              | Description                            |
|---------------------|----------------------------------------|
| sample_id           | Sample ID                              |
| bam_file            | Full path to `bam` file                |
| bam_bai_file        | Full path to `bam.bai` file            |

`--reference_genes_gtf_file`
* Reference genes GTF (`.gtf` or `.gtf.gz`) file.

`--params_longgf`
* Refer to the [LongGF documentation](https://github.com/WGLab/LongGF).

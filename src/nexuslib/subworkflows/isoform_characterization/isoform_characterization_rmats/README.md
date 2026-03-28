## isoform_characterization_rmats.nf

Characterizes isoforms using [rMATS](https://rnaseq-mats.sourceforge.io/).

### Inputs / Outputs

| I/O    | Description                                |
|:-------|:-------------------------------------------|
| Input  | `bam` and `bam.bai` files for each sample. | 
| Output | rMATS output files for each sample.        |

### Dependencies

* `rMATS`

### Example

```
nexus run --nf-workflow isoform_characterization_rmats.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --params_rmats '"-t paired --libType fr-unstranded --readLength 151"'
```

### Usage

```
workflow:
    1. Run rMATS.

usage: nexus run --nf-workflow isoform_characterization_rmats.nf [required] [optional] [--help]

required arguments:
    -c                             :   Nextflow .config file.
    -w                             :   Nextflow work directory path.
    --samples_tsv_file             :   TSV file with the following columns:
                                       'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                   :   Directory to which output files will be copied.
    --reference_genes_gtf_file     :   Reference genes GTF file.

optional arguments:
    --params_rmats                 :   ESPRESSO_S.pl parameters (default: '"-t paired --libType fr-unstranded --readLength 151"').
                                       Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                    |
|--------------|--------------------------------|
| sample_id    | Sample ID.                     |
| bam_file     | Full path to `bam` file.       |
| bam_bai_file | Full path to `bam.bai` file.   |

`--reference_genes_gtf_file`
* Reference genes GTF (`.gtf` or `.gtf.gz`) file.

`--params_rmats`
* Refer to the [rmats documentation](https://rnaseq-mats.sourceforge.io/).
* The following parameters for `rmats` are already included in `nexus` module for `rmats.py` and should not be specified:
  * `--b1`
  * `--gtf`
  * `--tmp`
  * `--nthread`
  * `--od`

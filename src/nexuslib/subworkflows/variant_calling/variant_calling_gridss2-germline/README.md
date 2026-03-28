## variant_calling_gridss.nf

Identifies somatic structural DNA variants in paired-read DNA BAM files using [GRIDSS](https://github.com/PapenfussLab/gridss).

### Inputs / Outputs

| I/O    | Description                                                |
|:-------|:-----------------------------------------------------------|
| Input  | Tumor and normal `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.                                |

### Dependencies

* `GRIDSS`

### Example

```
nexus run --nf-workflow variant_calling_gridss.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --params_gridss '""' \
    --params_gridss_somatic_filter '"--ref BSgenome.Hsapiens.UCSC.hg38"'
```

### Usage

```
workflow:
    1. Run GRIDSS (somatic mode).

usage: nexus run --nf-workflow variant_calling_gridss.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id',
                                            'tumor_bam_file',
                                            'tumor_bam_bai_file',
                                            'normal_bam_file',
                                            'normal_bam_bai_file'
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.

optional arguments:
    --params_gridss                     :   GRIDSS gridss parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
    --params_gridss_somatic_filter      :   GRIDSS gridss_somatic_filter parameters (default: '"--ref BSgenome.Hsapiens.UCSC.hg38"').
                                            Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header              | Description                        |
|---------------------|------------------------------------|
| sample_id           | Sample ID                          |
| tumor_bam_file      | Full path to tumor `bam` file      |
| tumor_bam_bai_file  | Full path to tumor `bam.bai` file  |
| normal_bam_file     | Full path to normal `bam` file     |
| normal_bam_bai_file | Full path to normal `bam.bai` file |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta`) file. An uncompressed `.fasta` file should be supplied.

`--params_gridss`
* Refer to the [GRIDSS documentation](https://github.com/PapenfussLab/gridss/blob/master/QuickStart.md).
* The following parameters for `gridss` are already included in `nexus` module and should not be specified:
  * `--reference`
  * `--output`
  * `--threads`
  * `--jvmheap`
  * `--otherjvmheap`

`--params_gridss_somatic_filter`
* Refer to the [GRIDSS documentation](https://github.com/PapenfussLab/gridss/blob/master/QuickStart.md).
* The following parameters for `gridss_somatic_filter` are already included in `nexus` module and should not be specified:
  * `--input`
  * `--output`
  * `--fulloutput`
  * `--scriptdir`
  * `-n`
  * `-t`

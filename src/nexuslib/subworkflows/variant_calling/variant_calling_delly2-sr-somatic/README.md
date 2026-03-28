## variant_calling_delly2.nf

Identifies somatic structural DNA variants in paired-end read DNA BAM files using [Delly2](https://github.com/dellytools/delly).

### Inputs / Outputs

| I/O    | Description                                                |
|:-------|:-----------------------------------------------------------|
| Input  | Tumor and normal `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.                                |

### Dependencies

* `Delly2`

### Example

```
nexus run --nf-workflow variant_calling_delly2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --exclude_tsv_file EXCLUDE_TSV_FILE \
    --params_delly2call '"--map-qual 20"'
```

### Usage

```
workflow:
    1. Run Delly2 (somatic mode).

usage: nexus run --nf-workflow variant_calling_delly2.nf [required] [optional] [--help]

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
    --exclude_tsv_file                  :   Delly2 'call' --exclude TSV file.

optional arguments:
    --params_delly2call                 :   Delly2 'call' parameters (default: '"--map-qual 20"').
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
* Reference genome FASTA (`.fasta` or `.fasta.gz`) file.

`--exclude_tsv_file`
* Delly2 'call' --exclude TSV file. This file can be downloaded from https://github.com/dellytools/delly/tree/main/excludeTemplates

`--params_delly2call`
* Refer to the [Delly2 documentation](https://github.com/dellytools/delly).
* The following parameters for `Delly2` are already included in `nexus` module for `delly call` and should not be specified:
  * `--outfile`
  * `--genome`
  * `--exclude`

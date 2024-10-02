## paired-end_read_dna_variant_calling_lumpy.nf

Identifies somatic structural DNA variants in paired-read DNA BAM files using [Lumpy](https://github.com/arq5x/lumpy-sv).

### Inputs / Outputs

| I/O    | Description                                                |
|:-------|:-----------------------------------------------------------|
| Input  | Tumor and normal `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.                                |

### Dependencies

* `Lumpy`

### Example

```
nexus run --nf-workflow paired-end_read_dna_variant_calling_lumpy.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR
```

### Usage

```
workflow:
    1. Run lumpyexpress (tumor and normal mode).

usage: nexus run --nf-workflow paired-end_read_dna_variant_calling_lumpy.nf [required] [optional] [--help]

required arguments:
    -c                                          :   Nextflow .config file.
    -w                                          :   Nextflow work directory path.
    --samples_tsv_file                          :   TSV file with the following columns:
                                                    'sample_id',
                                                    'tumor_bam_file',
                                                    'tumor_bam_bai_file',
                                                    'normal_bam_file',
                                                    'normal_bam_bai_file'
    --output_dir                                :   Directory to which output files will be copied.

optional arguments:
    --delete_work_dir                           :   Delete work directory (default: false).
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

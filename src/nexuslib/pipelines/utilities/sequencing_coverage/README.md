## sequencing_coverage.nf

Computes BAM sequencing coverage using [samtools](https://www.htslib.org/).

### Inputs / Outputs

| I/O    | Description                                |
|:-------|:-------------------------------------------|
| Input  | `bam` and `bam.bai` files for each sample. | 
| Output | `txt` file for each sample.                |

### Dependencies

* `samtools`

### Example

```
nexus run --nf-workflow sequencing_coverage.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR
```

### Usage

```
workflow:
    1. Run samtools coverage.

usage: nexus run --nf-workflow sequencing_coverage.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns: 
                                            'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --min_mapping_quality               :   Minimum mapping quality (default: 20).
    --min_base_quality                  :   Minimum base quality (default: 20).
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                     |
|--------------|---------------------------------|
| sample_id    | Sample ID.                      |
| bam_file     | Full path to `bam` file.        |
| bam_bai_file | Full path to `bam.bai` file.    |

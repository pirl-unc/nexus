## utilities_filter-rnabloom2-transcripts.nf

Filter RNA-Bloom2 output transcripts.

### Inputs / Outputs

| I/O    | Description                              |
|:-------|:-----------------------------------------|
| Input  | `rnabloom2_output_dir` each sample.      | 
| Output | `tsv` and `fasta` files for each sample. |

### Dependencies

* `nexus`

### Example

```
nexus run --nf-workflow utilities_filter-rnabloom2-transcripts.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --params_filter '"--min-mapping-quality 30 --min-read-support 3 --min-fraction-match 0.5"'
```

### Usage

```
workflow:
    1. Run filter_rnabloom2_transcripts (Nexus).

usage: nexus run --nf-workflow utilities_filter-rnabloom2-transcripts.nf [required] [optional] [--help]

required arguments:
    -c                      :   Nextflow .config file.
    -w                      :   Nextflow work directory path.
    --samples_tsv_file      :   TSV file with the following columns: 'sample_id', 'rnabloom2_output_dir''.
    --output_dir            :   Directory to which output files will be copied.

optional arguments:
    --params_filter         :   filter_rnabloom2_transcripts parameters (default: '"--min-mapping-quality 30 --min-read-support 3 --min-fraction-match 0.5"').
                                Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                                   |
|--------------|-----------------------------------------------|
| sample_id    | Sample ID.                                    |
| rnabloom2_output_dir     | Full path to RNA-bloom2 outputs directory.    |

`--params_filter`
* Refer to `filter_rnabloom2_transcripts -h` in `nexus`

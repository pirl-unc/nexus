## fastqc.nf

Runs a quality control check using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### Inputs / Outputs

| I/O    | Description                         |
|:-------|:------------------------------------|
| Input  | `fastq.gz` file(s) for each sample. | 
| Output | `txt` file for each sample.         |

### Dependencies

* `FastQC`

### Example

```
nexus run --nf-workflow fastqc.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --read_type 'single-end' OR 'paired-end'
```

### Usage

```
workflow:
    1. Run FastQC.

usage: nexus run --nf-workflow fastqc.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns: 
                                            'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --read_type                         :   read type (default: 'single-end'). 
                                            Either 'single-end' or 'paired-end'.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                                                            |
|--------------|------------------------------------------------------------------------|
| sample_id    | Sample ID.                                                             |
| fastq_file_1 | Full path to `fastq.gz` file.                                          |
| fastq_file_2 | Full path to `fastq.gz` file. (if single-end, leave this column empty) |

`--read_type`
* Either `single-end` or `paired-end`

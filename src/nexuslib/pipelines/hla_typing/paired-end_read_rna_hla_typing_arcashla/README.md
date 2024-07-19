## paired-end_read_rna_hla_typing_arcashla.nf

Performs HLA typing using paired-end (Illumina) RNA reads using [arcasHLA](https://github.com/RabadanLab/arcasHLA).

### Inputs / Outputs

| I/O    | Description                       |
|:-------|:----------------------------------|
| Input  | `R1 and R2 fastq.gz` files for each sample.       | 
| Output | `json` file for each sample.      |

### Dependencies

* `arcasHLA`

### Example

```
nexus run --nf-workflow paired-end_read_rna_hla_typing_arcashla.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR
```

### Usage

```
workflow:
    1. Profile HLA alleles using paired-end RNA sequencing BAM files using arcasHLA.

usage: nexus run --nf-workflow paired-end_read_rna_hla_typing_arcashla.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`

| Header       | Description                     |
| ------------ |---------------------------------|
| sample_id    | Sample ID.                      |
| fastq_file_1 | Full path to R1 `fastq.gz` file |
| fastq_file_2 | Full path to R2 `fastq.gz` file |


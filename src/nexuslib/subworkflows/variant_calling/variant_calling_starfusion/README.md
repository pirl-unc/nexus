## variant_calling_starfusion.nf

Identifies fusion genes in paired-read RNA fastq.gz files using [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion).

### Inputs / Outputs

| I/O    | Description                                            |
|:-------|:-------------------------------------------------------|
| Input  | RNA-seq `fastq.gz` files (R1 and R2) for each sample.   | 
| Output | `vcf` file for each sample.                            |

### Dependencies

* `STAR-Fusion`

### Example

```
nexus run --nf-workflow variant_calling_starfusion.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --genome_lib_dir GENOME_LIB_DIR
```

### Usage

```
workflow:
    1. Run STAR-Fusion.

usage: nexus run --nf-workflow variant_calling_starfusion.nf [required] [optional] [--help]

required arguments:
    -c                      :   Nextflow .config file.
    -w                      :   Nextflow work directory path.
    --samples_tsv_file      :   TSV file with the following columns:
                                'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir            :   Directory to which output files will be copied.
    --genome_lib_dir        :   Genome lib directory.

optional arguments:
    --params_starfusion     :   STAR-Fusion parameters (default: '""').
                                Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                             |
|--------------|-----------------------------------------|
| sample_id    | Sample ID                               |
| fastq_file_1 | Full path to RNA-seq `fastq.gz` R1 file |
| fastq_file_2 | Full path to RNA-seq `fastq.gz` R2 file |

`--genome_lib_dir`
* Genome library (CTAT resource library) directory. This can be downloaded [here](https://github.com/STAR-Fusion/STAR-Fusion/releases).

`--params_starfusion`
* Refer to the [STAR-Fusion documentation](https://github.com/STAR-Fusion/STAR-Fusion/wiki).
* The following parameters for `STAR-Fusion` are already included in `nexus` module for `starfusion` and should not be specified:
  * `--left_fq`
  * `--right_fq`
  * `--genome_lib_dir`
  * `--CPU`
  * `--output_dir`

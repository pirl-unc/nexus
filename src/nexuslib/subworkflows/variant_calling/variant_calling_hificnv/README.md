## variant_calling_hificnv.nf

Identifies copy number alterations in long-read DNA BAM files using [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV).

### Inputs / Outputs

| I/O    | Description                                         |
|:-------|:----------------------------------------------------|
| Input  | `bam`,  `bam.bai`, and `vcf` files for each sample. | 
| Output | `vcf` file for each sample.                         |

### Dependencies

* `Dysgu`

### Example

```
nexus run --nf-workflow variant_calling_hificnv.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --exclude_bed_file EXCLUDE_BED_FILE \
    --params_hificnv '""'
```

### Usage

```
workflow:
    1. Run hificnv.

usage: nexus run --nf-workflow variant_calling_hificnv.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'bam_file', 'bam_bai_file', 'vcf_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --exclude_bed_file                  :   Exclude BED file.

optional arguments:
    --params_hificnv                    :   hificnv parameters (default: '""').
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
| vcf_file            | Full path to small variants `vcf` file |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta`) file. An uncompressed `.fasta` file should be supplied.

`--exclude_bed_file`
* A BED file of genomic regions to exclude in the copy number alteration calculation.

`--params_hificnv`
* Refer to the [HiFiCNV documentation](https://github.com/PacificBiosciences/HiFiCNV).
* The following parameters for `hificnv` are already included in `nexus` module for `hificnv` and should not be specified:
  * `--ref`
  * `--bam`
  * `--maf`
  * `--exclude`
  * `--output-prefix`
  * `--threads`

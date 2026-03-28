## variant_calling_strelka2.nf

Identifies small somatic DNA variants (SNVs and INDELs) in paired-end DNA sequencing 
(BAM) files using [Strelka2](https://github.com/Illumina/strelka).

### Inputs / Outputs

| I/O    | Description                                                 |
|:-------|:------------------------------------------------------------|
| Input  | Tumor and normal `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.                                 |

### Dependencies

* `Strelka2`

### Example

```
nexus run --nf-workflow variant_calling_strelka2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --params_strelka2 '""'
```

### Usage

```
workflow:
    1.  Run Strelka2.

usage: nexus run --nf-workflow variant_calling_strelka2.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'
    --output_dir                        :   Directory to which output files will be symlinked.
    --reference_genome_fasta_file       :   Reference genome FASTA file.

optional arguments:
    --params_strelka2                   :   Strelka2 parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header              | Description                           |
|---------------------|---------------------------------------|
| sample_id           | Sample ID.                            |
| tumor_bam_file      | Full path to tumor `bam` file.        |
| tumor_bam_bai_file  | Full path to tumor `bam.bai` file.    |
| normal_bam_file     | Full path to normal `bam` file.       |
| normal_bam_bai_file | Full path to normal `bam.bai` file.   |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta` or `.fasta.gz`) file.

`--params_strelka2`
* Refer to the [Strelka2 documentation](https://github.com/Illumina/strelka).
* The following parameters for `Strelka2` are already included in `nexus` module for `Strelka2` and should not be specified:
  * `--normalBam`
  * `--tumorBam`
  * `--referenceFasta`
  * `--runDir`

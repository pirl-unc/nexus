## variant_calling_dysgu-lr.nf

Identifies somatic structural DNA variants in tumor and normal long-read DNA BAM files using [Dysgu](https://github.com/kcleal/dysgu).

### Inputs / Outputs

| I/O    | Description                                                   |
|:-------|:--------------------------------------------------------------|
| Input  | Tumor and normal `bam` and `bam.bai` files for each sample.   | 
| Output | `vcf` file for each sample.                                   |

### Dependencies

* `Dysgu`

### Example

```
nexus run --nf-workflow variant_calling_dysgu-lr.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --params_dysgu_run '"--mode pacbio --min-support 3 --min-size 30 --mq 20"' \
    --params_dysgu_filter '"--support-fraction 0.05 --min-mapq 20"'
```

### Usage

```
workflow:
    1. Run Dysgu.

usage: nexus run --nf-workflow variant_calling_dysgu-lr.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.

optional arguments:
    --params_dysgu_run                  :   Dysgu run parameters (default: '"--mode pacbio-revio --min-support 3 --min-size 30 --mq 20"').
                                            Note that the parameters need to be wrapped in quotes.
    --params_dysgu_filter               :   Dysgu filter parameters (default: '"--support-fraction 0.05 --min-mapq 20"').
                                            Note that the parameters need to be wrapped in quotes..
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

`--params_dysgu_run`
* Refer to the [Dysgu documentation](https://github.com/kcleal/dysgu).
* The following parameters for `dysgu run` are already included in `nexus` module for `dysgu` and should not be specified:
  * `--out-format`
  * `--procs`

`--params_dysgu_filter`
* Refer to the [Dysgu documentation](https://github.com/kcleal/dysgu).
* The following parameters for `dysgu filter` are already included in `nexus` module for `dysgu` and should not be specified:
  * `--normal-vcf`
  * `--procs`

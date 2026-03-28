## variant_calling_savana.nf

Identifies somatic structural DNA variants in tumor and normal long-read DNA BAM files using [Savana](https://github.com/cortes-ciriano-lab/savana).

### Inputs / Outputs

| I/O    | Description                                                |
|:-------|:-----------------------------------------------------------|
| Input  | Tumor and normal `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.                                |

### Dependencies

* `Savana`

### Example

```
nexus run --nf-workflow variant_calling_savana.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --custom_params_file CUSTOM_PARAMS_FILE \
    --contigs_txt_file CONTIGS_TXT_FILE \
    --params_savana_run '"--length 30 --mapq 20 --min_support 3"' \
    --params_savana_classify '""'
```

### Usage

```
workflow:
    1. Run Savana.

usage: nexus run --nf-workflow variant_calling_savana.nf [required] [optional] [--help]

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
    --custom_params_file                :   Savana classify --custom_params file.
    --contigs_txt_file                  :   TXT file with each contig name in a new line.

optional arguments:
    --params_savana_run                 :   Savana run parameters (default: '"--length 30 --mapq 20 --min_support 3"').
                                            Note that the parameters need to be wrapped in quotes.
    --params_savana_classify            :   Savana classify parameters (default: '""').
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

`--custom_params_file`
* Savana `classify` command `--custom_params` file.

`--contigs_txt_file`
* Savana `run` command `--contigs` file.

`--params_savana_run`
* Refer to the [savana documentation](https://github.com/cortes-ciriano-lab/savana).
* The following parameters for `savana run` are already included in `nexus` module for `savana` and should not be specified:
  * `-t`
  * `-n`
  * `--ref`
  * `--threads`
  * `--outdir`
  * `--sample`
  * `--contigs`

`--params_savana_classify`
* Refer to the [savana documentation](https://github.com/cortes-ciriano-lab/savana).
* The following parameters for `savana classify` are already included in `nexus` module for `savana` and should not be specified:
  * `--vcf`
  * `--output`
  * `--custom_params`

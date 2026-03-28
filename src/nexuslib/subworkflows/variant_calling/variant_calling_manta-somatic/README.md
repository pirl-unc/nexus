## variant_calling_manta.nf

Identifies somatic structural DNA variants in paired-read DNA BAM files using [Manta](https://github.com/Illumina/manta).

### Inputs / Outputs

| I/O    | Description                                                 |
|:-------|:------------------------------------------------------------|
| Input  | Tumor and normal `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.                                 |

### Dependencies

* `Manta`

### Example

```
nexus run --nf-workflow variant_calling_manta.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --params_manta_config '""' \
    --params_manta_run '""'
```

### Usage

```
workflow:
    1. Run manta (tumor and normal mode).

usage: nexus run --nf-workflow variant_calling_manta.nf [required] [optional] [--help]

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

optional arguments:
    --params_manta_config               :   Manta configManta.py parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
    --params_manta_run                  :   Manta runWorkflow.py parameters (default: '""').
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

`--params_manta_config`
* Refer to the [Manta documentation](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md).
* The following parameters for `Manta` `configManta.py` are already included in `nexus` module and should not be specified:
  * `--tumorBam`
  * `--normalBam`
  * `--referenceFasta`
  * `--runDir`

`--params_manta_run`
* Refer to the [Manta documentation](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md).
* The following parameters for `Manta` `runWorkflow.py` are already included in `nexus` module and should not be specified:
  * `-j`
  * `-g`

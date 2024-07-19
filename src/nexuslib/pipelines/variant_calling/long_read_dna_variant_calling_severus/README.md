## long_read_dna_variant_calling_severus.nf

Identifies somatic structural DNA variants in tumor and normal long-read DNA BAM files using [Severus](https://github.com/KolmogorovLab/Severus).

### Inputs / Outputs

| I/O    | Description                                   |
|:-------|:----------------------------------------------|
| Input  | Tumor and normal `bam` files for each sample. | 
| Output | `vcf` file for each sample.                   |

### Dependencies

* `Severus`

### Example

```
nexus run --nf-workflow long_read_dna_variant_calling_severus.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --params_severus '"--min-support 3 --min-sv-size 30 --min-mapq 20 --output-read-ids --bp-cluster-size 50"'
```

### Usage

```
workflow:
    1. Run Severus.

usage: nexus run --nf-workflow long_read_dna_variant_calling_severus.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --params_severus                :   Severus parameters (default: '"--min-support 3 --min-sv-size 30 --min-mapq 20 --output-read-ids --bp-cluster-size 50"').
                                        Note that the parameters need to be wrapped in quotes.
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`

| Header              | Description                        |
|---------------------|------------------------------------|
| sample_id           | Sample ID                          |
| tumor_bam_file      | Full path to tumor `bam` file      |
| tumor_bam_bai_file  | Full path to tumor `bam.bai` file  |
| normal_bam_file     | Full path to normal `bam` file     |
| normal_bam_bai_file | Full path to normal `bam.bai` file |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--params_severus`
* Refer to the [Severus documentation](https://github.com/KolmogorovLab/Severus).
* The following parameters for `Severus` are already included in `nexus` module for `severus` and should not be specified:
  * `--target-bam`
  * `--control-bam`
  * `--threads`
  * `--out-dir`

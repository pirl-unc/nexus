## long_read_dna_variant_calling_clairs.nf

Identifies somatic small DNA variants in tumor and normal long-read DNA BAM files using [ClairS](https://github.com/HKU-BAL/ClairS).

### Inputs / Outputs

| I/O    | Description                                   |
|:-------|:----------------------------------------------|
| Input  | Tumor and normal `bam` files for each sample. | 
| Output | `vcf` file for each sample.                   |

### Dependencies

* `ClairS`

### Example

```
nexus run --nf-workflow long_read_dna_variant_calling_clairs.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --params_clairs '"--platform hifi_revio --enable_indel_calling"'
```

### Usage

```
workflow:
    1. Run ClairS.

usage: nexus run --nf-workflow long_read_dna_variant_calling_clairs.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --params_clairs                     :   ClairS parameters (default: '"--platform hifi_revio --enable_indel_calling"').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
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

`--reference_genome_fasta_fai_file`
* Reference genome FASTA.FAI files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--params_clairs`
* Refer to the [ClairS documentation](https://github.com/HKU-BAL/ClairS).
* The following parameters for `run_clairs` are already included in `nexus` module for `clairs` and should not be specified:
  * `--tumor_bam`
  * `--normal_bam`
  * `--ref_fn`
  * `--output_dir`
  * `--threads`

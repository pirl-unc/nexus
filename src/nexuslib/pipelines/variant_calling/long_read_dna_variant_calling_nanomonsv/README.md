## long_read_dna_variant_calling_nanomonsv.nf

Identifies somatic structural DNA variants in tumor and normal long-read DNA BAM files using [Nanomonsv](https://github.com/friend1ws/nanomonsv).

### Inputs / Outputs

| I/O    | Description                                                   |
|:-------|:--------------------------------------------------------------|
| Input  | Tumor and normal `bam` and `bam.bai` files for each sample.   | 
| Output | `vcf` file for each sample.                                   |

### Dependencies

* `Nanomonsv`

### Example

```
nexus run --nf-workflow long_read_dna_variant_calling_nanomonsv.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --params_nanomonsv_parse '""' \
    --params_nanomonsv_get '"--min_tumor_variant_read_num 3 --min_tumor_VAF 0.05 --max_control_variant_read_num 0 --max_control_VAF 0.00 --min_indel_size 30 --max_panel_read_num 0 --median_mapQ_thres 20 --qv25"'
```

### Usage

```
workflow:
    1. Run Nanomonsv.

usage: nexus run --nf-workflow long_read_dna_variant_calling_nanomonsv.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --params_nanomonsv_parse            :   Nanomonsv parse parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
    --params_nanomonsv_get              :   Nanomonsv get parameters (default: '"--min_tumor_variant_read_num 3 --min_tumor_VAF 0.05 --max_control_variant_read_num 0 --max_control_VAF 0.00 --min_indel_size 30 --max_panel_read_num 0 --median_mapQ_thres 20 --qv25"').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
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
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--reference_genome_fasta_fai_file`
* Reference genome FASTA.FAI files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--params_nanomonsv_parse`
* Refer to the [Nanomonsv documentation](https://github.com/friend1ws/nanomonsv).
* The following parameters for `nanomonsv parse` are already included in `nexus` module for `nanomonsv` and should not be specified:
  * `--reference_fasta`

`--params_nanomonsv_get`
* Refer to the [Nanomonsv documentation](https://github.com/friend1ws/nanomonsv).
* The following parameters for `nanomonsv get` are already included in `nexus` module for `nanomonsv` and should not be specified:
  * `--control_prefix`
  * `--control_bam`
  * `--processes`

## paired-end_read_dna_variant_calling_gridss.nf

Identifies somatic structural DNA variants in paired-read DNA BAM files using [GRIDSS](https://github.com/PapenfussLab/gridss).

### Inputs / Outputs

| I/O    | Description                                   |
|:-------|:----------------------------------------------|
| Input  | Tumor and normal `bam` files for each sample. | 
| Output | `vcf` file for each sample.                   |

### Dependencies

* `GRIDSS`

### Example

```
nexus run --nf-workflow paired-end_read_dna_variant_calling_gridss.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --params_gridss '""' \
    --params_gridss_somatic_filter '"--ref BSgenome.Hsapiens.UCSC.hg38"'
```

### Usage

```
workflow:
    1. Run GRIDSS (somatic mode).

usage: nexus run --nf-workflow paired-end_read_dna_variant_calling_gridss.nf [required] [optional] [--help]

required arguments:
    -c                                          :   Nextflow .config file.
    -w                                          :   Nextflow work directory path.
    --samples_tsv_file                          :   TSV file with the following columns:
                                                    'sample_id',
                                                    'tumor_bam_file',
                                                    'tumor_bam_bai_file',
                                                    'normal_bam_file',
                                                    'normal_bam_bai_file'
    --output_dir                                :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file               :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file           :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --params_gridss                             :   GRIDSS gridss parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
    --params_gridss_somatic_filter              :   GRIDSS gridss_somatic_filter parameters (default: '"--ref BSgenome.Hsapiens.UCSC.hg38"').
                                                    Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                           :   Delete work directory (default: false).
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

`--params_gridss`
* Refer to the [GRIDSS documentation](https://github.com/PapenfussLab/gridss/blob/master/QuickStart.md).
* The following parameters for `gridss` are already included in `nexus` module and should not be specified:
  * `--reference`
  * `--output`
  * `--threads`
  * `--jvmheap`
  * `--otherjvmheap`

`--params_gridss_somatic_filter`
* Refer to the [GRIDSS documentation](https://github.com/PapenfussLab/gridss/blob/master/QuickStart.md).
* The following parameters for `gridss_somatic_filter` are already included in `nexus` module and should not be specified:
  * `--input`
  * `--output`
  * `--fulloutput`
  * `--scriptdir`
  * `-n`
  * `-t`

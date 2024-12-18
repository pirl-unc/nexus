## paired-end_read_dna_variant_calling_delly2.nf

Identifies somatic structural DNA variants in paired-end read DNA BAM files using [Delly2](https://github.com/dellytools/delly).

### Inputs / Outputs

| I/O    | Description                                                |
|:-------|:-----------------------------------------------------------|
| Input  | Tumor and normal `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.                                |

### Dependencies

* `Delly2`

### Example

```
nexus run --nf-workflow paired-end_read_dna_variant_calling_delly2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --exclude_tsv_file EXCLUDE_TSV_FILE \
    --params_delly2call '"--map-qual 20"'
```

### Usage

```
workflow:
    1. Run Delly2 (somatic mode).

usage: nexus run --nf-workflow paired-end_read_dna_variant_calling_delly2.nf [required] [optional] [--help]

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
    --exclude_tsv_file                          :   Delly2 'call' --exclude TSV file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/delly2/human.hg38.excl.tsv).
    --params_delly2call                         :   Delly2 'call' parameters (default: '"--map-qual 20"').
                                                    Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                           :   Delete work directory (default: false).
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

`--exclude_tsv_file`
* Delly2 'call' --exclude TSV files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/delly2/ on LBG.

`--params_delly2call`
* Refer to the [Delly2 documentation](https://github.com/dellytools/delly).
* The following parameters for `Delly2` are already included in `nexus` module for `delly call` and should not be specified:
  * `--outfile`
  * `--genome`
  * `--exclude`
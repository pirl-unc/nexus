## long_read_dna_variant_calling_cutesv.nf

Identifies structural DNA variants in long-read DNA BAM files using [cuteSV](https://github.com/tjiangHIT/cuteSV).

### Inputs / Outputs

| I/O    | Description                                |
|:-------|:-------------------------------------------|
| Input  | `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.                |

### Dependencies

* `cuteSV`

### Example

```
nexus run --nf-workflow long_read_dna_variant_calling_cutesv.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --params_cutesv '"--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_support 3 --min_mapq 20 --min_size 30 --max_size -1 --report_readid --genotype"'
```

### Usage

```
workflow:
    1. Run CuteSV.

usage: nexus run --nf-workflow long_read_dna_variant_calling_cutesv.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --params_cutesv                 :   CuteSV parameters (default:
                                        '"--max_cluster_bias_INS 1000
                                          --diff_ratio_merging_INS 0.9
                                          --max_cluster_bias_DEL 1000
                                          --diff_ratio_merging_DEL 0.5
                                          --min_support 3
                                          --min_mapq 20
                                          --min_size 30
                                          --max_size -1
                                          --report_readid
                                          --genotype"').
                                        Note that the parameters need to be wrapped in quotes.
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                 |
|--------------|-----------------------------|
| sample_id    | Sample ID                   |
| bam_file     | Full path to `bam` file     |
| bam_bai_file | Full path to `bam.bai` file |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--params_cutesv`
* Refer to the [cuteSV documentation](https://github.com/tjiangHIT/cuteSV).
* The following parameters for `cuteSV` are already included in `nexus` module for `cuteSV` and should not be specified:
  * `--threads`
  * `--sample`

## long_read_rna_variant_calling_clair3rna.nf

Identifies RNA variants in long-read RNA BAM files using [Clair3-RNA](https://github.com/HKU-BAL/Clair3-RNA).

### Inputs / Outputs

| I/O    | Description                                |
|:-------|:-------------------------------------------|
| Input  | `bam` and `bam.bai` files for each sample. | 
| Output | Clair3-RNA outputs                         |

### Dependencies

* `Clair3-RNA`

### Example

```
nexus run --nf-workflow long_read_rna_variant_calling_clair3rna.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --params_clair3rna '"--platform hifi_sequel2_minimap2"'
```

### Usage

```
workflow:
    1. Run Clair3-RNA.

usage: nexus run --nf-workflow long_read_rna_variant_calling_clair3rna.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --params_clair3rna                  :   Clair3-RNA parameters (default: '"--platform hifi_sequel2_minimap2"').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                  |
|--------------|------------------------------|
| sample_id    | Sample ID.                   |
| bam_file     | Full path to `bam` file.     |
| bam_bai_file | Full path to `bam.bai` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--reference_genome_fasta_fai_file`
* Reference genome FASTA.FAI files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--params_clair3rna`
* Refer to the [Clair3-RNA documentation](https://github.com/HKU-BAL/Clair3-RNA).
* The following parameters for `Clair3-RNA` are already included in `nexus` module for `run_clair3_rna` and should not be specified:
  * `--bam_fn`
  * `--ref_fn`
  * `--output_dir`
  * `--threads`

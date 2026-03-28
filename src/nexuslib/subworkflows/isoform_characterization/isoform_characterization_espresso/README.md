## isoform_characterization_espresso.nf

Characterizes isoforms using [espresso](https://github.com/Xinglab/espresso).

### Inputs / Outputs

| I/O    | Description                                 |
|:-------|:--------------------------------------------|
| Input  | `bam` and `bam.bai` files for each sample.  | 
| Output | Espresso outputs for each sample.           |

### Dependencies

* `espresso`

### Example

```
nexus run --nf-workflow isoform_characterization_espresso.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --params_espresso_s '""' \
    --params_espresso_c '""' \
    --params_espresso_q '""'
```

### Usage

```
workflow:
    1. Run ESPRESSO.

usage: nexus run --nf-workflow isoform_characterization_espresso.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file (uncompressed .fasta file).
    --reference_genes_gtf_file          :   Reference genes GTF file (uncompressed .gtf file).

optional arguments:
    --params_espresso_s                 :   ESPRESSO_S.pl parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
    --params_espresso_c                 :   ESPRESSO_C.pl parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
    --params_espresso_q                 :   ESPRESSO_Q.pl parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header     | Description                   |
| ---------- |-------------------------------|
| sample_id  | Sample ID.                    |
| fastq_file | Full path to `fastq.gz` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta`) file. An uncompressed `.fasta` file should be supplied.

`--reference_genes_gtf_file`
* Reference genes GTF (`.gtf`) file. An uncompressed `.gtf` file should be supplied.

`--params_espresso_s`
* Refer to the [espresso documentation](https://github.com/Xinglab/espresso).
* The following parameters for `ESPRESSO_S.pl` are already included in `nexus` module for `espresso` and should not be specified:
  * `--list_samples`
  * `--fa`
  * `--anno`
  * `--out`
  * `--num_thread`

`--params_espresso_c`
* Refer to the [espresso documentation](https://github.com/Xinglab/espresso).
* The following parameters for `ESPRESSO_C.pl` are already included in `nexus` module for `espresso` and should not be specified:
  * `--in`
  * `--fa`
  * `--target_ID`
  * `--num_thread`

`--params_espresso_q`
* Refer to the [espresso documentation](https://github.com/Xinglab/espresso).
* The following parameters for `ESPRESSO_Q.pl` are already included in `nexus` module for `espresso` and should not be specified:
  * `--list_samples`
  * `--anno`
  * `--num_thread`

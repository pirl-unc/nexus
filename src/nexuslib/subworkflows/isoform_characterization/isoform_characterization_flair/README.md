## isoform_characterization_flair.nf

Characterizes isoforms using [flair](https://github.com/BrooksLabUCSC/flair).

### Inputs / Outputs

| I/O    | Description                          |
|:-------|:-------------------------------------|
| Input  | `fastq.gz` file for each sample.     | 
| Output | Flair output files for each sample.  |

### Dependencies

* `flair`

### Example

```
nexus run --nf-workflow isoform_characterization_flair.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --params_flair_align '""' \
    --params_flair_correct '""' \
    --params_flair_collapse '""'
```

### Usage

```
workflow:
    1. Run flair align.
    2. Run flair correct.
    3. Run flair collapse.

usage: nexus run --nf-workflow isoform_characterization_flair.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --reference_genes_gtf_file          :   Reference genes GTF file.

optional arguments:
    --params_flair_align                :   Flair 'align' parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
    --params_flair_correct              :   Flair 'correct' parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
    --params_flair_collapse             :   Flair 'collapse' parameters (default: '" "').
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

`--params_flair_align`
* Refer to the [flair documentation](https://github.com/BrooksLabUCSC/flair).
* The following parameters for `flair` are already included in `nexus` module for `flair align` and should not be specified:
  * `-g`
  * `-r`
  * `--output`
  * `--threads`

`--params_flair_correct`
* Refer to the [flair documentation](https://github.com/BrooksLabUCSC/flair).
* The following parameters for `flair` are already included in `nexus` module for `flair correct` and should not be specified:
  * `--query`
  * `--genome`
  * `--gtf`
  * `--output`
  * `--threads`

`--params_flair_collapse`
* Refer to the [flair documentation](https://github.com/BrooksLabUCSC/flair).
* The following parameters for `flair` are already included in `nexus` module for `flair collapse` and should not be specified:
  * `--query`
  * `--genome`
  * `--gtf`
  * `--reads`
  * `--output`
  * `--threads`

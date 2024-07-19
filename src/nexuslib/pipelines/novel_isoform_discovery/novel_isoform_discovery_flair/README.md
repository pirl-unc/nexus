## novel_isoform_discovery_flair.nf

Discovers novel isoforms using [flair](https://github.com/BrooksLabUCSC/flair).

### Inputs / Outputs

| I/O    | Description                                                                                 |
|:-------|:--------------------------------------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                                                            | 
| Output | Flair outputs including `fasta`, `gtf`, `bed`, `bam`, and `bam.bai` files  for each sample. |

### Dependencies

* `flair`

### Example

```
nexus run --nf-workflow novel_isoform_discovery_flair.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --gtf_file GTF_FILE
```

### Usage

```
workflow:
    1. Run flair align.
    2. Run flair correct.
    3. Run flair collapse.

usage: nexus run --nf-workflow novel_isoform_discovery_flair.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --gtf_file                      :   Reference transcriptome GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf).
    --params_flair_align            :   flair 'align' parameters (default: '""').
                                        Note that the parameters need to be wrapped in quotes.
    --params_flair_correct          :   flair 'correct' parameters (default: '""').
                                        Note that the parameters need to be wrapped in quotes.
    --params_flair_collapse         :   flair 'collapse' parameters (default: '""').
                                        Note that the parameters need to be wrapped in quotes.
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`

| Header     | Description                   |
| ---------- |-------------------------------|
| sample_id  | Sample ID.                    |
| fastq_file | Full path to `fastq.gz` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

`--gtf_file`
* Reference GTF files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

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

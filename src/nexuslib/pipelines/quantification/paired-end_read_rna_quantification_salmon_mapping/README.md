## paired-end_read_rna_quantification_salmon_mapping.nf

Quantify RNA abundance in paired-end RNA reads using [Salmon](https://salmon.readthedocs.io/en/latest/index.html) (mapping mode).

### Inputs / Outputs

| I/O    | Description                                                                              |
|:-------|:-----------------------------------------------------------------------------------------|
| Input  | `R1 and R2 fastq.gz` files for each sample.                                              | 
| Output | Salmon output files including `quant.sf` and `quant.genes.sf` files for each sample. |

### Dependencies

* `Salmon`

### Example

```
nexus run --nf-workflow paired-end_read_rna_quantification_salmon_mapping.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_transcripts_fasta_file REFERENCE_TRANSCRIPTS_FASTA_FILE \
    --params_salmon_index '"--gencode"' \
    --params_salmon_quant '"--libType IU --seqBias --gcBias --posBias"'
```

### Usage

```
workflow:
    1. Quantify RNA in paired-end FASTQ files using Salmon.

usage: nexus run --nf-workflow paired-end_read_rna_quantification_salmon_mapping.nf [required] [optional] [--help]

required arguments:
    -c                                      :   Nextflow .config file.
    -w                                      :   Nextflow work directory path.
    --samples_tsv_file                      :   TSV file with the following columns:
                                                'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir                            :   Directory to which output files will be copied.

optional arguments:
    --reference_transcripts_fasta_file      :   Reference transcripts FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.transcripts.fa).
    --params_salmon_index                   :   Salmon index parameters (default: '"--gencode"').
                                                Note that the parameters need to be wrapped in quotes.
    --params_salmon_quant                   :   Salmon parameters (default: '"--libType IU --seqBias --gcBias --posBias"').
                                                Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                       :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                     |
| ------------ |---------------------------------|
| sample_id    | Sample ID.                      |
| fastq_file_1 | Full path to R1 `fastq.gz` file |
| fastq_file_2 | Full path to R2 `fastq.gz` file |

`--reference_transcripts_fasta_file`
* Reference transcripts FASTA files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

`--params_salmon_index`
* Refer to the [Salmon documentation](https://salmon.readthedocs.io/en/latest/index.html).
* The following parameters for `Salmon` are already included in `nexus` module for `salmon index` and should not be specified:
  * `--transcripts`
  * `--index`
  * `--threads`

`--params_salmon_quant`
* Refer to the [Salmon documentation](https://salmon.readthedocs.io/en/latest/index.html).
* The following parameters for `Salmon` are already included in `nexus` module for `salmon quant` and should not be specified:
  * `--index`
  * `--mates1`
  * `--mates2`
  * `--geneMap`
  * `--output`
  * `--threads`

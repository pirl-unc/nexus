## long_read_rna_quantification_oarfish.nf

Quantify RNA abundance in long RNA reads using [Oarfish](https://github.com/COMBINE-lab/oarfish).

### Inputs / Outputs

| I/O    | Description                                                   |
|:-------|:--------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                              | 
| Output | Oarfish output files including `.quant` file for each sample. |

### Dependencies

* `Oarfish`

### Example

```
nexus run --nf-workflow long_read_rna_quantification_oarfish.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --reference_transcriptome_fasta_file REFERENCE_TRANSCRIPTOME_FASTA_FILE \
    --output_dir OUTPUT_DIR \
    --params_oarfish '"--seq-tech pac-bio-hifi"'
```

### Usage

```
workflow:
    1. Quantify RNA in long-read FASTQ files using Oarfish.

usage: nexus run --nf-workflow long_read_rna_quantification_oarfish.nf [required] [optional] [--help]

required arguments:
    -c                                      :   Nextflow .config file.
    -w                                      :   Nextflow work directory path.
    --samples_tsv_file                      :   TSV file with the following columns:
                                                'sample_id', 'fastq_file'.
    --reference_transcriptome_fasta_file    :   Reference transcriptome FASTA file
                                                (default: '/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.transcripts.fa').
    --params_oarfish                        :   Oarfish parameters (default: '"--seq-tech pac-bio-hifi"').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
    --output_dir                            :   Directory to which output files will be copied.

optional arguments:
    --delete_work_dir                       :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                   |
|--------------|-------------------------------|
| sample_id    | Sample ID.                    |
| fastq_file   | Full path to `fastq.gz` file. |

`--reference_transcriptome_fasta_file`
* Reference transcriptome FASTA files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

`--params_oarfish`
* Refer to the [Oarfish documentation](https://github.com/COMBINE-lab/oarfish).
* The following parameters for `Oarfish` are already included in `nexus` module for `oarfish` and should not be specified:
  * `--reference`
  * `--threads`
  * `--output`
  * `--reads`


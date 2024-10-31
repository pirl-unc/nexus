## long_read_rna_quantification_transigner.nf

Quantify RNA abundance in long RNA reads using [TranSigner](https://github.com/haydenji0731/transigner).

### Inputs / Outputs

| I/O    | Description                                                                                     |
|:-------|:------------------------------------------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                                                                | 
| Output | TranSigner output files including `abundances.tsv` and `assignments.tsv` files for each sample. |

### Dependencies

* `TranSigner`

### Example

```
nexus run --nf-workflow long_read_rna_quantification_transigner.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --reference_transcriptome_fasta_file REFERENCE_TRANSCRIPTOME_FASTA_FILE \
    --output_dir OUTPUT_DIR \
    --params_transigner_align '""' \
    --params_transigner_prefilter '"--filter -tp -500 -fp -600"' \
    --params_transigner_em '"--drop --use-score"'
```

### Usage

```
workflow:
    1. Quantify RNA in long-read FASTQ files using TranSigner.

usage: nexus run --nf-workflow long_read_rna_quantification_transigner.nf [required] [optional] [--help]

required arguments:
    -c                                      :   Nextflow .config file.
    -w                                      :   Nextflow work directory path.
    --samples_tsv_file                      :   TSV file with the following columns:
                                                'sample_id', 'fastq_file'.
    --reference_transcriptome_fasta_file    :   Reference transcriptome FASTA file
                                                (default: '/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.transcripts.fa').
    --params_transigner_align               :   TranSigner align parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
    --params_transigner_prefilter           :   TranSigner prefilter parameters (default: '"--filter -tp -500 -fp -600"').
                                                Note that the parameters need to be wrapped in quotes
                                                and a space at the end of the string is necessary.
    --params_transigner_em                  :   TranSigner em parameters (default: '"--drop --use-score"').
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

`--params_transigner_align`
* Refer to the [TranSigner documentation](https://github.com/haydenji0731/transigner).
* The following parameters for `transigner align` are already included in `nexus` module for `transigner` and should not be specified:
  * `-q`
  * `-t`
  * `-d`
  * `-o`
  * `-p`

`--params_transigner_prefilter`
* Refer to the [TranSigner documentation](https://github.com/haydenji0731/transigner).
* The following parameters for `transigner prefilter` are already included in `nexus` module for `transigner` and should not be specified:
  * `-a`
  * `-t`
  * `-o`

`--params_transigner_em`
* Refer to the [TranSigner documentation](https://github.com/haydenji0731/transigner).
* The following parameters for `transigner em` are already included in `nexus` module for `transigner` and should not be specified:
  * `-s`
  * `-i`
  * `-o`

  


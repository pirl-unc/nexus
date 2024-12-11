## peptide_alignment_blastp.nf

Align peptide sequences to a reference proteome using [blastp](https://blast.ncbi.nlm.nih.gov/doc/blast-help/).

### Inputs / Outputs

| I/O    | Description                                   |
|:-------|:----------------------------------------------|
| Input  | `fasta.gz` file for each sample.              | 
| Output | Blastp `-outfmt 7` `txt` file for each sample.  |

### Dependencies

* `blast`

### Example
```
nexus run --nf-workflow peptide_alignment_minimap2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_proteome_fasta_file REFERENCE_PROTEOME_FASTA_FILE \
    --params_blastp '""'
```

### Usage

```
workflow:
    1. Align peptide sequences (fasta.gz files) to a reference proteome using Blastp.

usage: nexus run --nf-workflow peptide_alignment_blastp.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fasta_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_proteome_fasta_file     :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.pc_translations.fa).
    --params_blastp                     :   Blastp parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow `config` file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header     | Description                   |
|------------|-------------------------------|
| sample_id  | Sample ID.                    |
| fasta_file | Full path to `fasta.gz` file. |

`--reference_proteome_fasta_file`
* Reference proteome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--params_blastp`
* Refer to the [blastp documentation](https://blast.ncbi.nlm.nih.gov/doc/blast-help/).
* The following parameters for `blastp` are already included in `nexus` module for `blastp` and should not be specified:
  * `-query`
  * `-db`
  * `-out`
  * `-outfmt`
  * `-num_threads`

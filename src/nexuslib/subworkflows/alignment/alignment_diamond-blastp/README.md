## alignment_diamond-blastp.nf

Run Blastp using [DIAMOND](https://github.com/bbuchfink/diamond).

### Inputs / Outputs

| I/O    | Description                                |
|:-------|:-------------------------------------------|
| Input  | `fasta.gz` file for each sample.           | 
| Output | DIAMOND Blastp `tsv` file for each sample. |

### Dependencies

* `diamond`

### Example
```
nexus run --nf-workflow alignment_diamond-blastp.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --database_dmnd_file DATABASE_DNMD_FILE \
    --params_diamond_blastp '""'
```

### Usage

```
workflow:
    1. Blast peptide sequences (fasta or fasta.gz files) using Diamond Blastp.

usage: nexus run --nf-workflow alignment_diamond-blastp.nf [required] [optional] [--help]

required arguments:
    -c                          :   Nextflow .config file.
    -w                          :   Nextflow work directory path.
    --samples_tsv_file          :   TSV file with the following columns:
                                    'sample_id', 'fasta_file'.
    --output_dir                :   Directory to which output files will be copied.
    --database_dmnd_file        :   Diamond-prepared local Blast database file.

optional arguments:
    --params_diamond_blastp     :   Diamond blastp parameters (default: '""').
                                    Note that the parameters need to be wrapped in quotes
                                    and a space at the end of the string is necessary.
```

### Parameters

`-c`
* Nextflow `config` file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header     | Description                              |
|------------|------------------------------------------|
| sample_id  | Sample ID.                               |
| fasta_file | Full path to `fasta` or `fasta.gz` file. |

`--database_dmnd_file`
* Diamond-prepared local Blast database file.
* To prepare a DMND file, first download ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz. Then run `diamond makedb --threads 64 --in nr.gz -d nr`.

`--params_diamond_blastp`
* Refer to the [blastp documentation](https://blast.ncbi.nlm.nih.gov/doc/blast-help/).
* The following parameters for `diamond blastp` are already included in `nexus` module for `diamond blastp` and should not be specified:
  * `--query`
  * `--db`
  * `--out`
  * `--header`
  * `--verbose`
  * `--threads`

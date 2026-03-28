## alignment_blastp.nf

Query peptide sequences using [blastp](https://blast.ncbi.nlm.nih.gov/doc/blast-help/).

### Inputs / Outputs

| I/O    | Description                                   |
|:-------|:----------------------------------------------|
| Input  | `fasta.gz` file for each sample.              | 
| Output | Blastp `-outfmt 7` `txt` file for each sample.  |

### Dependencies

* `blast`

### Example
```
nexus run --nf-workflow alignment_blastp.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --blastdb_dir BLASTDB_DIR \
    --blastdb_name BLASTDB_NAME \
    --params_blastp '""'
```

### Usage

```
workflow:
    1. Blast peptide sequences (fasta or fasta.gz files) using Blastp.

usage: nexus run --nf-workflow alignment_blastp.nf [required] [optional] [--help]

required arguments:
    -c                      :   Nextflow .config file.
    -w                      :   Nextflow work directory path.
    --samples_tsv_file      :   TSV file with the following columns:
                                'sample_id', 'fasta_file'.
    --output_dir            :   Directory to which output files will be copied.
    --blastdb_dir           :   Local Blast database path. The database can be downloaded by the
                                following command: update_blastdb.pl --source gcp --decompress nr --num_threads 16

optional arguments:
    --blastdb_name          :   Local Blast database name (default: 'nr').
    --params_blastp         :   Blastp parameters (default: '""').
                                Note that the parameters need to be wrapped in quotes
                                and a space at the end of the string is necessary.
```

### Parameters

`-c`
* Nextflow `config` file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header     | Description                               |
|------------|-------------------------------------------|
| sample_id  | Sample ID.                                |
| fasta_file | Full path to `fastq` or `fasta.gz` file.  |

`--blastdb_dir`
* The database can be downloaded by: `update_blastdb.pl --source gcp --decompress nr --num_threads 16`

`--blastdb_name`
* Local Blast database name (`nr` by default).

`--params_blastp`
* Refer to the [blastp documentation](https://blast.ncbi.nlm.nih.gov/doc/blast-help/).
* The following parameters for `blastp` are already included in `nexus` module for `blastp` and should not be specified:
  * `-query`
  * `-db`
  * `-out`
  * `-outfmt`
  * `-num_threads`

## peptide_alignment_blastp.nf

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
nexus run --nf-workflow peptide_alignment_diamond_blastp.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --blastdb_dir BLASTDB_DIR \
    --blastdb_name BLASTDB_NAME \
    --params_diamond_blastp '""'
```

### Usage

```
workflow:
    1. Run blastp using DIAMOND.

usage: nexus run --nf-workflow peptide_alignment_diamond_blastp.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fasta_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --blastdb_dir                       :   Diamond-prepared local Blast database path (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/blast/).
    --blastdb_name                      :   Local Blast database name (default: 'nr').
    --params_diamond_blastp             :   Diamond blastp parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
    --delete_work_dir                   :   Delete work directory (default: false).
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

`--blastdb_dir`
* Diamond-prepared local Blast database path (`/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/blast/` on LBG).
* The database can be downloaded and prepared by: 
  * `update_blastdb.pl --source gcp --decompress nr --num_threads 16`
  * `diamond prepdb --db nr --threads 16`

`--blastdb_name`
* Local Blast database name (`nr` on LBG).
* The database can be downloaded and prepared by: 
  * `update_blastdb.pl --source gcp --decompress nr --num_threads 16`
  * `diamond prepdb --db nr --threads 16`

`--params_diamond_blastp`
* Refer to the [blastp documentation](https://blast.ncbi.nlm.nih.gov/doc/blast-help/).
* The following parameters for `diamond blastp` are already included in `nexus` module for `diamond blastp` and should not be specified:
  * `--query`
  * `--db`
  * `--out`
  * `--header`
  * `--verbose`
  * `--threads`

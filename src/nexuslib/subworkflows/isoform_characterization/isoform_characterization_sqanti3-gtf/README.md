## isoform_characterization_sqanti3-gtf.nf

Characterizes isoforms using [Sqanti3](https://github.com/ConesaLab/SQANTI3).

### Inputs / Outputs

| I/O    | Description                           |
|:-------|:--------------------------------------|
| Input  | `gtf` file for each sample.           | 
| Output | Sqanti3 output files for each sample. |

### Dependencies

* `Sqanti3`

### Example

```
nexus run --nf-workflow isoform_characterization_sqanti3-gtf.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --params_sqanti3_qc '"--report skip"' \
    --params_sqanti3_filter '"ml"'
```

### Usage

```
workflow:
    1. Run SQANTI3 (GTF mode).

usage: nexus run --nf-workflow isoform_characterization_sqanti3-gtf.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'gtf_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --reference_genes_gtf_file          :   Reference genes GTF file.

optional arguments:
    --params_sqanti3_qc                 :   sqanti3_qc.py parameters (default: '"--report skip"').
                                            Note that the parameters need to be wrapped in quotes.
    --params_sqanti3_filter             :   sqanti3_filter.py parameters (default: '"ml"').
                                            Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header    | Description              |
|-----------|--------------------------|
| sample_id | Sample ID.               |
| gtf_file  | Full path to `gtf` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta`) file. An uncompressed `.fasta` file should be supplied.

`--reference_genes_gtf_file`
* Reference genes GTF (`.gtf`) file. An uncompressed `.gtf` file should be supplied.

`--params_sqanti3_qc`
* Refer to the [Sqanti3 documentation](https://github.com/ConesaLab/SQANTI3).
* The following parameters for `sqanti3` are already included in `nexus` module for `sqanti3_qc.py` and should not be specified:
  * `--fasta`
  * `-t`
  * `-o`
  * `-d`

`--params_sqanti3_filter`
* Refer to the [Sqanti3 documentation](https://github.com/ConesaLab/SQANTI3).
* The following parameters for `sqanti3` are already included in `nexus` module for `sqanti3_filter.py` and should not be specified:
  * `--output`
  * `--dir`

## isoform_characterization_mandalorion.nf

Characterizes isoforms using [Mandalorion](https://github.com/christopher-vollmers/Mandalorion).

### Inputs / Outputs

| I/O    | Description                               |
|:-------|:------------------------------------------|
| Input  | `fastq.gz` file for each sample.          | 
| Output | Mandalorion output files for each sample. |

### Dependencies

* `Mandalorion`

### Example

```
nexus run --nf-workflow isoform_characterization_mandalorion.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --params_mandalorion '""'
```

### Usage

```
workflow:
    1. Run Mandalorion.

usage: nexus run --nf-workflow isoform_characterization_mandalorion.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --reference_genes_gtf_file          :   Reference genes GTF file.

optional arguments:
    --params_mandalorion                :   Mando.py parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header     | Description                   |
|------------|-------------------------------|
| sample_id  | Sample ID.                    |
| fastq_file | Full path to `fastq.gz` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta` or `.fasta.gz`) file.

`--reference_genes_gtf_file`
* Reference genes GTF (`.gtf` or `.gtf.gz`) file.

`--params_mandalorion`
* Refer to the [mandalorion documentation](https://github.com/christopher-vollmers/Mandalorion).
* The following parameters for `mandalorion` are already included in `nexus` module for `Mando.py` and should not be specified:
  * `--path`
  * `--genome_annotation`
  * `--genome_sequence`
  * `-f`
  * `-t`

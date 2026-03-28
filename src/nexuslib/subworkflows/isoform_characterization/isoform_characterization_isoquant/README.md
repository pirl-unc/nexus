## isoform_characterization_isoquant.nf

Characterizes isoforms using [IsoQuant](https://github.com/ablab/IsoQuant).

### Inputs / Outputs

| I/O    | Description                            |
|:-------|:---------------------------------------|
| Input  | `fastq.gz` file for each sample.       | 
| Output | IsoQuant output files for each sample. |

### Dependencies

* `IsoQuant`

### Example

```
nexus run --nf-workflow isoform_characterization_isoquant.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --params_isoquant '"--data_type pacbio_ccs --sqanti_output --high_memory --complete_genedb"'
```

### Usage

```
workflow:
    1. Run isoquant.

usage: nexus run --nf-workflow isoform_characterization_isoquant.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --reference_genes_gtf_file          :   Reference genes GTF file.

optional arguments:
    --params_isoquant                   :   isoquant parameters (default: '"--data_type pacbio_ccs --sqanti_output --high_memory --complete_genedb"').
                                            Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header     | Description                   |
| ---------- |-------------------------------|
| sample_id  | Sample ID.                    |
| fastq_file | Full path to `fastq.gz` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta`) file. An uncompressed `.fasta` file should be supplied.

`--reference_genes_gtf_file`
* Reference genes GTF (`.gtf`) file. An uncompressed `.gtf` file should be supplied.

`--params_isoquant`
* Refer to the [IsoQuant documentation](https://github.com/ablab/IsoQuant).
* The following parameters for `IsoQuant` are already included in `nexus` module for `IsoQuant` and should not be specified:
  * `--reference`
  * `--genedb`
  * `--fastq_list`
  * `--threads`
  * `-o`

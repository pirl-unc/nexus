## quantification_bambu.nf

Quantify RNA abundance in long RNA reads using [Bambu](https://github.com/GoekeLab/bambu).

### Inputs / Outputs

| I/O    | Description                                |
|:-------|:-------------------------------------------|
| Input  | `bam` and `bam.bai` files for each sample. | 
| Output | Bambu output files                         |

### Dependencies

* `Bambu`

### Example

```
nexus run --nf-workflow quantification_bambu.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --output_dir OUTPUT_DIR
```

### Usage

```
workflow:
    1. Quantify RNA in long-read FASTQ files using Bambu.

usage: nexus run --nf-workflow quantification_bambu.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                    :   Directory to which output files will be copied.
    --reference_genome_fasta_file   :   Reference genome FASTA file.
    --reference_genes_gtf_file      :   Reference genes GTF file.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                  |
|--------------|------------------------------|
| sample_id    | Sample ID.                   |
| bam_file     | Full path to `bam` file.     |
| bam_bai_file | Full path to `bam.bai` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta` or `fasta.gz`) file.

`--reference_genes_gtf_file`
* Reference genes GTF (`.gtf` or `.gtf.gz`) file.

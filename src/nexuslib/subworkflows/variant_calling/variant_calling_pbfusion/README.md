## variant_calling_pbfusion.nf

Identifies fusion genes in long-read RNA BAM files using [pbfusion](https://github.com/PacificBiosciences/pbfusion).

### Inputs / Outputs

| I/O    | Description                               |
|:-------|:------------------------------------------|
| Input  | `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.               |

### Dependencies

* `pbfusion`

### Example

```
nexus run --nf-workflow variant_calling_pbfusion.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genes_gtf_file REFERENCE_GTF_FILE \
    --params_pbfusion_discover '"--min-coverage 3 --min-mean-mapq 20 --gtf-transcript-allow-lncRNA"'
```

### Usage

```
workflow:
    1. Run Pbfusion.

usage: nexus run --nf-workflow variant_calling_pbfusion.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id',
                                            'bam_file',
                                            'bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --reference_genes_gtf_file          :   Reference genes GTF file.

optional arguments:
    --params_pbfusion_discover          :   Pbfusion discover parameters (default: '"--min-coverage 3 --min-mean-mapq 20 --gtf-transcript-allow-lncRNA"').
                                            Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header              | Description                            |
|---------------------|----------------------------------------|
| sample_id           | Sample ID                              |
| bam_file            | Full path to `bam` file                |
| bam_bai_file        | Full path to `bam.bai` file            |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta` or `.fasta.gz`) file.

`--reference_genes_gtf_file`
* Reference genes GTF (`.gtf` or `.gtf.gz`) file.

`--params_pbfusion_discover`
* Refer to the [pbfusion documentation](https://github.com/PacificBiosciences/pbfusion).

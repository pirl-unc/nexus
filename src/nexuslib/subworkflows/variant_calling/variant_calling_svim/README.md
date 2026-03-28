## variant_calling_svim.nf

Identifies structural DNA variants in long-read DNA BAM files using [SVIM](https://github.com/eldariont/svim).

### Inputs / Outputs

| I/O    | Description                                |
|:-------|:-------------------------------------------|
| Input  | `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.                |

### Dependencies

* `SVIM`

### Example

```
nexus run --nf-workflow variant_calling_svim.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --params_svim '"--min_mapq 20 --min_sv_size 30 --insertion_sequences --read_names --zmws"'
```

### Usage

```
workflow:
    1. Run SVIM.

usage: nexus run --nf-workflow variant_calling_svim.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.

optional arguments:
    --params_svim                       :   SVIM parameters (default: '"--min_mapq 20 --min_sv_size 30 --insertion_sequences --read_names --zmws"').
                                            Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                 |
|--------------|-----------------------------|
| sample_id    | Sample ID                   |
| bam_file     | Full path to `bam` file     |
| bam_bai_file | Full path to `bam.bai` file |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta` or `.fasta.gz`) file.

`--params_svim`
* Refer to the [SVIM documentation](https://github.com/eldariont/svim).
* The following parameters for `SVIM` are already included in `nexus` module for `svim` and should not be specified:
  * `--sample`

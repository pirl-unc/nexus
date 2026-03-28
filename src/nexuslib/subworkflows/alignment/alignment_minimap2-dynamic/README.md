## alignment_minimap2-dynamic.nf

Aligns long DNA or RNA reads (Pacific Biosciences or Oxford Nanopore) to a reference genome using [minimap2](https://github.com/lh3/minimap2).

### Inputs / Outputs

| I/O    | Description                                                     |
|:-------|:----------------------------------------------------------------|
| Input  | `fastq.gz` file and `fasta` file for each sample.               | 
| Output | MD-tagged and sorted `bam` and `bam.bai` files for each sample. |

### Dependencies

* `minimap2`
* `samtools`

### Example
```
nexus run --nf-workflow alignment_minimap2-dynamic.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --params_minimap2 '"-ax map-hifi --cs --eqx -Y -L --secondary=no"'
```

### Usage

```
workflow:
    1. Align reads (fastq.gz files) to a custom reference genome using minimap2.
    2. Convert sam files to bam files.
    3. Generate MD tags.
    4. Sort MD-tagged bam files.

usage: nexus run --nf-workflow alignment_minimap2-dynamic.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file', 'reference_genome_fasta_file', 'reference_genome_fasta_fai_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --params_minimap2                   :   Minimap2 parameters (default: "-ax map-hifi --cs --eqx -Y -L --secondary=no").
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
    --platform_tag                      :   Platform tag (default: 'unknown').
    --platform_unit_tag                 :   Platform unit tag (default: 'unknown').
    --library_tag                       :   Library tag (default: 'unknown').
```

### Parameters

`-c`
* Nextflow `config` file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header                          | Description                    |
|---------------------------------|--------------------------------|
| sample_id                       | Sample ID.                     |
| fastq_file                      | Full path to `fastq.gz` file.  |
| reference_genome_fasta_file     | Full path to `fasta` file.     |
| reference_genome_fasta_fai_file | Full path to `fasta.fai` file. |

`--params_minimap2`
* Refer to the [minimap2 documentation](https://github.com/lh3/minimap2).
* The following parameters for `minimap2` are already included in `nexus` module for `minimap2` and should not be specified:
  * `-t`
  * `-R`

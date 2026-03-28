## alignment_ultra.nf

Aligns long RNA reads (Pacific Biosciences or Oxford Nanopore) to a reference genome using [uLTRA](https://github.com/ksahlin/ultra).

### Inputs / Outputs

| I/O    | Description                                                     |
|:-------|:----------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                                | 
| Output | MD-tagged and sorted `bam` and `bam.bai` files for each sample. |

### Dependencies

* `uLTRA`
* `samtools`

### Example

```
nexus run --nf-workflow alignment_ultra.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --params_ultra_index '"--disable_infer"' \
    --params_ultra '"--isoseq"'
```

### Usage

```
workflow:
    1. Align reads (fastq.gz files) to a reference genome using uLTRA.
    2. Generate MD tags.
    3. Sort MD-tagged bam file.

usage: nexus run --nf-workflow alignment_ultra.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --reference_genes_gtf_file          :   Reference genes GTF file.

optional arguments:
    --params_ultra_index                :   uLTRA index parameters (default: "--disable_infer").
    --params_ultra                      :   uLTRA align parameters (default: "--isoseq").
```

### Parameters

`-c`
* Nextflow `config` file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header     | Description                  |
| ---------- |------------------------------|
| sample_id  | Sample ID.                   |
| fastq_file | Full path to `fastq.gz` file |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta` or `.fasta.gz`) file.

`--reference_genes_gtf_file`
* Reference genes GTF (`.gtf` or `.gtf.gz`) file.

`--params_ultra_index`
* `uLTRA index` parameters.
* Refer to the [uLTRA documentation](https://github.com/ksahlin/ultra).

`--params_ultra`
* `uLTRA align` parameters.
* Refer to the [uLTRA documentation](https://github.com/ksahlin/ultra).
* The following parameters for `minimap2` are already included in `nexus` module for `minimap2` and should not be specified:
  * `--t`
  * `--index`
  * `--prefix`

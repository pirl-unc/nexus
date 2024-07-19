## long_read_alignment_minimap2.nf

Aligns long DNA or RNA reads (Pacific Biosciences or Oxford Nanopore) to a reference genome using [minimap2](https://github.com/lh3/minimap2).

### Inputs / Outputs

| I/O    | Description                                                      |
|:-------|:-----------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                                 | 
| Output | MD-tagged and sorted `bam` and `bam.bai` files for each sample.  |

### Dependencies

* `minimap2`
* `samtools`

### Example
```
nexus run --nf-workflow long_read_alignment_minimap2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --params_minimap2 '"-ax map-hifi --cs --eqx -Y -L"'
```

### Usage

```
workflow:
    1. Align reads (fastq.gz files) to a reference genome using minimap2.
    2. Convert sam files to bam files.
    3. Generate MD tags.
    4. Sort MD-tagged bam files.

usage: nexus run --nf-workflow long_read_alignment_minimap2.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --params_minimap2                   :   Minimap2 parameters (default: "-ax map-hifi --cs --eqx -Y -L").
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
    --platform_tag                      :   Platform tag (default: 'unknown').
    --platform_unit_tag                 :   Platform unit tag (default: 'unknown').
    --library_tag                       :   Library tag (default: 'unknown').
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow `config` file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header     | Description                   |
| ---------- |-------------------------------|
| sample_id  | Sample ID.                    |
| fastq_file | Full path to `fastq.gz` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--reference_genome_fasta_fai_file`
* Reference genome FASTA.FAI files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--params_minimap2`
* Refer to the [minimap2 documentation](https://github.com/lh3/minimap2).
* The following parameters for `minimap2` are already included in `nexus` module for `minimap2` and should not be specified:
  * `-t`
  * `-R`